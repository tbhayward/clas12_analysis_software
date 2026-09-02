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
#include <TPad.h>

#include <nlohmann/json.hpp>

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cctype>
#include <cstdio>
#include <cstdlib>
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
#include <unordered_set>
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


static void configure_production_pad(TPad* pad,
                                     double left = 0.14,
                                     double right = 0.035,
                                     double bottom = 0.14,
                                     double top = 0.075) {
    if (!pad) return;
    pad->SetLeftMargin(left);
    pad->SetRightMargin(right);
    pad->SetBottomMargin(bottom);
    pad->SetTopMargin(top);
    pad->SetTicks(1, 1);
    pad->SetGrid(0, 0);
}

static void style_production_frame(TH1* h,
                                   double x_title_size = 0.050,
                                   double y_title_size = 0.050,
                                   double label_size = 0.040) {
    if (!h) return;
    h->SetStats(0);
    h->SetLineColor(kBlack);
    h->SetLineWidth(1);
    h->GetXaxis()->SetTitleFont(42);
    h->GetYaxis()->SetTitleFont(42);
    h->GetXaxis()->SetLabelFont(42);
    h->GetYaxis()->SetLabelFont(42);
    h->GetXaxis()->SetTitleSize(x_title_size);
    h->GetYaxis()->SetTitleSize(y_title_size);
    h->GetXaxis()->SetLabelSize(label_size);
    h->GetYaxis()->SetLabelSize(label_size);
    h->GetXaxis()->SetTitleOffset(1.10);
    h->GetYaxis()->SetTitleOffset(1.28);
}

static std::unique_ptr<TLatex> make_period_label(const std::string& period,
                                                 double x = 0.94,
                                                 double y = 0.88,
                                                 int align = 33,
                                                 double size = 0.050) {
    auto label = std::make_unique<TLatex>();
    label->SetNDC();
    label->SetTextAlign(align);
    label->SetTextFont(62);
    label->SetTextSize(size);
    label->DrawLatex(x, y, period.c_str());
    return label;
}


static std::pair<double, double> production_efficiency_plot_range(
    const std::string& quantity) {
    // Production plotting convention: a pathological low-statistics point may
    // extend beyond the visible frame, but it must never determine the scale of
    // the entire canvas.  Numerical values remain fully available in CSV.
    if (quantity == "data_over_mc" || quantity == "relative_ratio") {
        return {0.50, 1.50};
    }
    return {0.20, 1.60};
}

static double production_pull_limit() {
    // Keep closure canvases readable.  Pulls outside this range remain in CSV.
    return 6.0;
}

static bool is_polar_angle_model_variable(const std::string& key) {
    return key == "e_theta" || key == "p_theta" || key == "g_theta";
}

static std::string pretty_kinematic_key(const std::string& key) {
    if (key == "Q2") return "Q^{2}";
    if (key == "xB") return "x_{B}";
    if (key == "t") return "|t|";
    if (key == "phi") return "#phi";
    if (key == "e_theta") return "#theta_{e}";
    if (key == "p_theta") return "#theta_{p}";
    if (key == "g_theta") return "#theta_{#gamma}";
    return key;
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

    try {
        fin >> j;
    } catch (const std::exception& e) {
        fatal("[current_dependence] FATAL: failed to parse combined cuts JSON " +
              path + ": " + e.what());
    }

    if (!j.is_object()) {
        fatal("[current_dependence] FATAL: combined cuts JSON is not an object: " +
              path);
    }

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
            static const std::set<std::string> supported_variables = {
                "Delta_phi",
                "theta",
                "theta_gamma_gamma",
                "pTmiss",
                "Emiss2",
                "Mx2",
                "Mx2_2"
            };

            for (const auto& kv : vm) {
                if (supported_variables.find(kv.first) ==
                    supported_variables.end()) {
                    fatal("[current_dependence] FATAL: cut key '" + key +
                          "' contains unsupported production exclusivity "
                          "variable '" + kv.first + "'.");
                }
            }

            // The Python optimizer always starts with Mx2. Other variables are
            // optional because the greedy optimization stops when no further
            // cut improves S/sqrt(S+B).
            if (vm.find("Mx2") == vm.end()) {
                fatal("[current_dependence] FATAL: cut key '" + key +
                      "' is missing the mandatory first production cut 'Mx2'.");
            }

            out.emplace(key, std::move(vm));
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
        ss << "[current_dependence] FATAL: production exclusivity cut requires "
           << "a missing ROOT branch for logical variable '" << var << "'.";
        if (var == "theta_gamma_gamma") {
            ss << " Direct ep->eppi0 samples require physical branch "
               << "'theta_pi0_pi0'.";
        }
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
static std::mutex g_progress_output_mutex;

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

    double theta = 0.0; bool has_theta = false;
    double theta_gamma_gamma = 0.0; bool has_theta_gamma_gamma = false;
    double theta_pi0_pi0 = 0.0; bool has_theta_pi0_pi0 = false;
    double Emiss2 = 0.0; bool has_Emiss2 = false;
    double Mx2 = 0.0; bool has_Mx2 = false;
    double Mx2_2 = 0.0; bool has_Mx2_2 = false;

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
        // Python-optimized production exclusivity variables. The direct-pi0
        // JSON uses the logical theta_gamma_gamma variable, while its physical
        // ROOT branch is theta_pi0_pi0.
        ena("Delta_phi");
        ena("theta");
        ena("theta_gamma_gamma");
        ena("theta_pi0_pi0");
        ena("pTmiss");
        ena("Emiss2");
        ena("Mx2");
        ena("Mx2_2");

        ena("e_p");
        ena("e_theta");
        ena("e_phi");

        ena("p1_theta");
        ena("p1_phi");

        ena("p2_p");
        ena("p2_theta");
        ena("p2_phi");

        t->SetCacheSize(16LL * 1024LL * 1024LL);
        t->AddBranchToCache("*", true);

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

        bD("theta", &theta, has_theta);
        bD("theta_gamma_gamma", &theta_gamma_gamma, has_theta_gamma_gamma);
        bD("theta_pi0_pi0", &theta_pi0_pi0, has_theta_pi0_pi0);
        bD("Emiss2", &Emiss2, has_Emiss2);
        bD("Mx2", &Mx2, has_Mx2);
        bD("Mx2_2", &Mx2_2, has_Mx2_2);

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

    if (global_cuts_require_auxiliary_kinematics(cfg)) {
        if (!(b.has_e_theta && b.has_e_phi &&
              b.has_p1_theta && b.has_p1_phi &&
              b.has_p2_p && b.has_p2_theta && b.has_p2_phi)) {
            fatal("[current_dependence] FATAL: auxiliary fiducial cuts require e_theta, e_phi, p1_theta, p1_phi, p2_p, p2_theta, p2_phi branches.");
        }

        return passes_global_cuts(b.t1, b.open_angle_ep2, b.pTmiss,
                                  b.detector1, b.detector2,
                                  tags.period_label,
                                  b.e_p, b.e_theta, b.e_phi,
                                  b.p1_theta, b.p1_phi,
                                  b.p2_p, b.p2_theta, b.p2_phi,
                                  cfg);
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
                                      b.e_p, b.e_theta, b.e_phi, b.p1_theta, b.p1_phi,
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
                                  tags.period_label,
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

    bool has_delta_phi = false;
    const double delta_phi = b.delta_phi_value(has_delta_phi);

    if (!check_sigma_var(vm, "Delta_phi",
                         has_delta_phi, delta_phi)) {
        return false;
    }

    if (!check_sigma_var(vm, "theta",
                         b.has_theta, b.theta)) {
        return false;
    }

    // The Python output keeps theta_gamma_gamma as the logical variable name
    // for both channels. Direct ep->eppi0 trees store the corresponding angle
    // in theta_pi0_pi0.
    if (cfg.channel == Channel::EPPI0) {
        if (!check_sigma_var(vm, "theta_gamma_gamma",
                             b.has_theta_pi0_pi0,
                             b.theta_pi0_pi0)) {
            return false;
        }
    } else {
        if (!check_sigma_var(vm, "theta_gamma_gamma",
                             b.has_theta_gamma_gamma,
                             b.theta_gamma_gamma)) {
            return false;
        }
    }

    if (!check_sigma_var(vm, "pTmiss",
                         b.has_pTmiss, b.pTmiss)) {
        return false;
    }

    if (!check_sigma_var(vm, "Emiss2",
                         b.has_Emiss2, b.Emiss2)) {
        return false;
    }

    if (!check_sigma_var(vm, "Mx2",
                         b.has_Mx2, b.Mx2)) {
        return false;
    }

    if (!check_sigma_var(vm, "Mx2_2",
                         b.has_Mx2_2, b.Mx2_2)) {
        return false;
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

struct ResolvedRunInfo {
    bool valid = false;
    int current = 0;
    double charge = std::numeric_limits<double>::quiet_NaN();
};

static ResolvedRunInfo resolve_run_info_cached(
    const PeriodTags& tags,
    int runnum,
    const std::unordered_map<int, ChargeEntry>& charge_map,
    bool use_second_column_charge_for_all_unpolarized,
    bool use_columns_3_to_5_charge_sum_scaled_for_fa18_sp19_unpolarized,
    double columns_3_to_5_charge_sum_scale) {

    ResolvedRunInfo out;

    if (!resolve_current(tags.internal, runnum, out.current)) {
        return out;
    }

    const auto charge_it = charge_map.find(runnum);
    if (charge_it == charge_map.end()) {
        return out;
    }

    if (!select_unpolarized_charge_for_period(
            tags, runnum, charge_it->second,
            use_second_column_charge_for_all_unpolarized,
            use_columns_3_to_5_charge_sum_scaled_for_fa18_sp19_unpolarized,
            columns_3_to_5_charge_sum_scale, out.charge)) {
        return out;
    }

    out.valid = out.charge > 0.0;
    return out;
}

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
    std::unordered_map<int, ResolvedRunInfo> run_info_cache;
    run_info_cache.reserve(256);

    for (Long64_t i = 0; i < N; ++i) {
        tree->GetEntry(i);

        if (!passes_cone_cut(b)) {
            continue;
        }

        auto cache_it = run_info_cache.find(b.runnum);
        if (cache_it == run_info_cache.end()) {
            cache_it = run_info_cache.emplace(
                b.runnum,
                resolve_run_info_cached(
                    tags, b.runnum, charge_map,
                    use_second_column_charge_for_all_unpolarized,
                    use_columns_3_to_5_charge_sum_scaled_for_fa18_sp19_unpolarized,
                    columns_3_to_5_charge_sum_scale)).first;
        }

        const ResolvedRunInfo& run_info = cache_it->second;
        if (!run_info.valid) {
            continue;
        }

        if (!passes_global_dispatch(b, tags)) {
            continue;
        }

        if (!passes_sigma_dispatch(cfg, tags, data_cuts, b)) {
            continue;
        }

        out.counts_by_run[b.runnum] += 1;
        out.current_by_run[b.runnum] = run_info.current;
        out.charge_by_run[b.runnum] = run_info.charge;
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

static double charge_weighted_data_current(const std::vector<CurrentPoint>& points) {
    // The scan ordinate is yield/charge, so the response appropriate to an
    // integrated production sample must be averaged with accumulated charge,
    // not with the already efficiency-biased reconstructed event counts.
    double total_charge = 0.0;
    double weighted_current = 0.0;

    for (const CurrentPoint& p : points) {
        if (std::isfinite(p.charge) && p.charge > 0.0) {
            total_charge += p.charge;
            weighted_current += p.charge * (double)p.current_nA;
        }
    }

    if (!(total_charge > 0.0)) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    return weighted_current / total_charge;
}

static double weighted_data_rel(const std::vector<CurrentPoint>& points,
                                const FitResult& fit) {
    const double weighted_current = charge_weighted_data_current(points);
    if (!std::isfinite(weighted_current)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return rel_at_current(weighted_current, fit);
}

static double weighted_data_rel_err(const std::vector<CurrentPoint>& points,
                                    const FitResult& fit) {
    const double weighted_current = charge_weighted_data_current(points);
    if (!std::isfinite(weighted_current)) {
        return std::numeric_limits<double>::quiet_NaN();
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
    lat.DrawLatex(0.18, 0.50, "Sp19 Inb efficiency");
    lat.DrawLatex(0.18, 0.42, "copied from Fa18 Inb");
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
        lat.DrawLatex(0.16, 0.76, "Sp19 Inb efficiency copied from Fa18 Inb");
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
        lat.DrawLatex(0.15, 0.61, "Sp19 Inb efficiency copied from Fa18 Inb");
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
                                      bool use_fa18_for_sp19,
                                      bool draw_legend = true) {
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
            lat.DrawLatex(0.08, 0.17, "Sp19 Inb efficiency copied from Fa18 Inb");
        }
    }

    if (draw_legend) {
        leg_top->Draw();
    }

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


static void draw_data_mc_period_legend_panel(
    const std::vector<PeriodResult>& results,
    bool use_fa18_for_sp19) {

    configure_production_pad((TPad*)gPad, 0.08, 0.05, 0.08, 0.08);

    TH1F* frame = gPad->DrawFrame(0.0, 0.0, 1.0, 1.0);
    frame->SetStats(0);
    frame->GetXaxis()->SetLabelSize(0.0);
    frame->GetYaxis()->SetLabelSize(0.0);
    frame->GetXaxis()->SetTickLength(0.0);
    frame->GetYaxis()->SetTickLength(0.0);

    TLatex lat;
    lat.SetNDC(true);
    lat.SetTextFont(62);
    lat.SetTextSize(0.052);
    lat.DrawLatex(0.10, 0.88, "Run periods");

    double y = 0.74;
    for (const std::string& period : PERIOD_ORDER) {
        if (period == "Sp19 Inb" && use_fa18_for_sp19) continue;
        auto it = std::find_if(results.begin(), results.end(),
                               [&](const PeriodResult& r) { return r.period == period; });
        if (it == results.end()) continue;

        lat.SetTextColor(period_color(period));
        lat.SetTextFont(62);
        lat.SetTextSize(0.041);
        lat.DrawLatex(0.15, y, period.c_str());
        y -= 0.10;
    }

    lat.SetTextColor(kBlack);
    lat.SetTextFont(42);
    lat.SetTextSize(0.031);
    lat.DrawLatex(0.10, 0.26, "filled markers: DATA");
    lat.DrawLatex(0.10, 0.20, "open markers: MC");

    if (use_fa18_for_sp19) {
        lat.SetTextColor(kRed + 2);
        lat.SetTextFont(42);
        lat.SetTextSize(0.029);
        lat.DrawLatex(0.10, 0.10, "Sp19 Inb efficiency copied from Fa18 Inb");
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
            // The unused Sp19 panel is the cleanest place for the shared
            // run-period legend and the concise Sp19 transfer note.
            draw_data_mc_period_legend_panel(results, use_fa18_for_sp19);
        } else if (it != results.end()) {
            draw_data_mc_period_stack(
                &(*it), period, false, results, use_fa18_for_sp19, true);
        } else {
            draw_data_mc_period_stack(
                nullptr, period, false, results, use_fa18_for_sp19, true);
        }

        ++pad;
    }

    c.cd(6);
    draw_data_mc_period_stack(
        nullptr, "All periods", true, results, use_fa18_for_sp19, false);

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
        lat.DrawLatex(0.15, 0.84, "Sp19 Inb efficiency copied from Fa18 Inb.");
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
        lat.DrawLatex(0.15, 0.84, "Sp19 Inb efficiency copied from Fa18 Inb.");
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
            lat.DrawLatex(0.16, 0.58, "Sp19 Inb efficiency copied from Fa18 Inb");
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
    // Sp19 lacks a usable DATA luminosity scan, so transfer only the DATA
    // current response from Fa18 Inb.  Keep the independently determined Sp19
    // MC response: its current-overlay reconstruction study is valid and need
    // not inherit the DATA limitation.  Carry the DATA fit itself because the
    // event-level correction needs its relative slope m/b.
    dst.data_factor = src.data_factor;
    dst.data_factor_err = src.data_factor_err;
    dst.data_fit = src.data_fit;

    std::cout << "[current_dependence] " << reason
              << ": copied DATA current response from "
              << src.period << " to " << dst.period
              << " data_factor=" << dst.data_factor
              << " +/- " << dst.data_factor_err
              << "; retained " << dst.period
              << " MC factor=" << dst.mc_factor
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
    std::unordered_map<int, ResolvedRunInfo> run_info_cache;
    run_info_cache.reserve(256);

    for (Long64_t i = 0; i < N; ++i) {
        tree->GetEntry(i);

        if (!passes_cone_cut(b)) {
            continue;
        }

        auto cache_it = run_info_cache.find(b.runnum);
        bool new_run = false;
        if (cache_it == run_info_cache.end()) {
            cache_it = run_info_cache.emplace(
                b.runnum,
                resolve_run_info_cached(
                    tags, b.runnum, charge_map,
                    use_second_column_charge_for_all_unpolarized,
                    use_columns_3_to_5_charge_sum_scaled_for_fa18_sp19_unpolarized,
                    columns_3_to_5_charge_sum_scale)).first;
            new_run = true;
        }

        const ResolvedRunInfo& run_info = cache_it->second;
        if (!run_info.valid) {
            continue;
        }

        if (new_run) {
            for (const KinematicVarConfig& v : vars) {
                for (DataAgg& agg : out[v.key]) {
                    agg.counts_by_run[b.runnum] += 0;
                    agg.current_by_run[b.runnum] = run_info.current;
                    agg.charge_by_run[b.runnum] = run_info.charge;
                }
            }
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
            agg.current_by_run[b.runnum] = run_info.current;
            agg.charge_by_run[b.runnum] = run_info.charge;
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


static TGraphErrors* make_region_theta_factor_graph(const std::vector<KinematicBinResult>& bins,
                                                    int color) {
    TGraphErrors* g = new TGraphErrors();
    int ip = 0;

    for (const KinematicBinResult& b : bins) {
        // Diagnostic plotting policy: do not impose an arbitrary minimum number
        // of current settings, event count, relative uncertainty, or slope sign.
        // Two-current-point fits are therefore shown.  The only plotting sanity
        // requirement is that the fitted efficiency itself be finite and
        // physically positive with a finite non-negative uncertainty.  Every
        // fitted bin, including pathological/non-positive fits, remains in CSV.
        const bool usable =
            std::isfinite(b.factor) &&
            std::isfinite(b.factor_err) &&
            b.factor > 0.0 &&
            b.factor_err >= 0.0;
        if (!usable) continue;

        g->SetPoint(ip, b.x_center, b.factor);
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


static std::vector<KinematicBinResult> make_data_over_mc_kinematic_bins(
    const std::vector<KinematicBinResult>& data_bins,
    const std::vector<KinematicBinResult>& mc_bins) {

    std::vector<KinematicBinResult> out;
    const size_t n = std::min(data_bins.size(), mc_bins.size());
    out.reserve(n);

    for (size_t i = 0; i < n; ++i) {
        const KinematicBinResult& d = data_bins[i];
        const KinematicBinResult& m = mc_bins[i];

        KinematicBinResult r;
        r.x_low = d.x_low;
        r.x_high = d.x_high;
        r.x_center = std::isfinite(d.x_center) ? d.x_center : m.x_center;
        r.x_err = d.x_err;
        r.n_current_points = std::min(d.n_current_points, m.n_current_points);
        r.total_counts = std::min(d.total_counts, m.total_counts);

        if (std::isfinite(d.factor) && d.factor > 0.0 &&
            std::isfinite(m.factor) && m.factor > 0.0 &&
            std::isfinite(d.factor_err) && d.factor_err >= 0.0 &&
            std::isfinite(m.factor_err) && m.factor_err >= 0.0) {

            r.factor = d.factor / m.factor;
            const double rel2 =
                std::pow(d.factor_err / d.factor, 2) +
                std::pow(m.factor_err / m.factor, 2);
            r.factor_err = std::fabs(r.factor) * std::sqrt(std::max(0.0, rel2));
        }

        out.push_back(r);
    }

    return out;
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

    auto draw_info_panel =
        [&](const std::string& sample_label, bool fit_version) {
        configure_production_pad((TPad*)gPad, 0.08, 0.05, 0.08, 0.08);
        TH1F* info_frame = gPad->DrawFrame(0.0, 0.0, 1.0, 1.0);
        info_frame->SetStats(0);
        info_frame->GetXaxis()->SetLabelSize(0.0);
        info_frame->GetYaxis()->SetLabelSize(0.0);
        info_frame->GetXaxis()->SetTickLength(0.0);
        info_frame->GetYaxis()->SetTickLength(0.0);

        TLatex lat;
        lat.SetNDC(true);
        lat.SetTextFont(62);
        lat.SetTextSize(0.050);
        lat.DrawLatex(0.10, 0.89,
                      (cfg.title + " " + sample_label).c_str());

        lat.SetTextFont(42);
        lat.SetTextSize(0.034);
        lat.DrawLatex(0.10, 0.81, "Current-efficiency dependence");

        double y = 0.69;
        for (const std::string& period : PERIOD_ORDER) {
            if (hide_sp19_inb_from_replacement_plots(
                    hide_sp19_inb_from_all_period_plots, period)) continue;
            lat.SetTextColor(period_color(period));
            lat.SetTextFont(62);
            lat.SetTextSize(0.035);
            lat.DrawLatex(0.13, y, period.c_str());
            y -= 0.075;
        }
        lat.SetTextColor(kBlack);
        lat.SetTextFont(42);
        lat.SetTextSize(0.029);
        if (fit_version) {
            lat.DrawLatex(0.10, 0.28, "Curves: weighted linear fits");
            lat.DrawLatex(0.10, 0.22, "Bands: 1#sigma fit uncertainty");
        } else {
            lat.DrawLatex(0.10, 0.25, "Points: current-efficiency factors");
        }
        lat.DrawLatex(0.10, 0.14, "Same event selection in all panels");
    };

    auto draw_one_canvas = [&](const std::string& canvas_name,
                               const std::string& output_name,
                               const std::string& sample_label,
                               const VarPeriodResults& results_by_var_period,
                               const std::map<std::string, std::map<std::string, LinearFitSummary>>& fits_by_var_period,
                               const std::map<std::string, double>& integrated,
                               bool fit_version,
                               const std::string& y_label) {
        (void)integrated;
        TCanvas c(canvas_name.c_str(), canvas_name.c_str(), 2400, 1200);
        c.Divide(4, 2, 0.002, 0.002);

        std::vector<std::unique_ptr<TObject>> owned;

        for (size_t iv = 0; iv < vars.size(); ++iv) {
            const KinematicVarConfig& v = vars[iv];
            c.cd((int)iv + 1);
            configure_production_pad((TPad*)gPad, 0.16, 0.035, 0.16, 0.08);

            const double xmin = v.edges.front();
            const double xmax = v.edges.back();

            const auto yrange = production_efficiency_plot_range(
                sample_label == "DATA/MC" ? "data_over_mc" : "factor");
            const double ymin = yrange.first;
            const double ymax = yrange.second;

            auto frame = std::make_unique<TH1D>(
                ("h_kinematic_prod_" + canvas_name + "_" + std::to_string(iv)).c_str(),
                (";" + v.x_label + ";" + y_label).c_str(),
                100, xmin, xmax);
            frame->SetDirectory(nullptr);
            frame->SetMinimum(ymin);
            frame->SetMaximum(ymax);
            style_production_frame(frame.get(), 0.055, 0.055, 0.043);
            frame->Draw();
            owned.emplace_back(std::move(frame));

            auto unity = std::make_unique<TLine>(xmin, 1.0, xmax, 1.0);
            unity->SetLineStyle(2);
            unity->SetLineWidth(1);
            unity->SetLineColor(kGray + 1);
            unity->Draw("SAME");
            owned.emplace_back(std::move(unity));

            // Draw shaded fit bands first, then fit lines, then points.
            if (fit_version) {
                for (const std::string& period : PERIOD_ORDER) {
                    if (hide_sp19_inb_from_replacement_plots(
                            hide_sp19_inb_from_all_period_plots, period)) continue;
                    const auto it_var = fits_by_var_period.find(v.key);
                    if (it_var == fits_by_var_period.end()) continue;
                    const auto it_fit = it_var->second.find(period);
                    if (it_fit == it_var->second.end() || !it_fit->second.valid) continue;
                    TGraph* band =
                        make_linear_fit_band_graph(
                            it_fit->second, xmin, xmax, period_color(period));
                    band->Draw("F SAME");
                }

                for (const std::string& period : PERIOD_ORDER) {
                    if (hide_sp19_inb_from_replacement_plots(
                            hide_sp19_inb_from_all_period_plots, period)) continue;
                    const auto it_var = fits_by_var_period.find(v.key);
                    if (it_var == fits_by_var_period.end()) continue;
                    const auto it_fit = it_var->second.find(period);
                    if (it_fit == it_var->second.end() || !it_fit->second.valid) continue;
                    TGraph* line =
                        make_linear_fit_line_graph(
                            it_fit->second, xmin, xmax, period_color(period));
                    line->SetLineWidth(2);
                    line->Draw("L SAME");
                }
            }

            for (const std::string& period : PERIOD_ORDER) {
                if (hide_sp19_inb_from_replacement_plots(
                        hide_sp19_inb_from_all_period_plots, period)) continue;
                if (it_var_results == results_by_var_period.end()) continue;
                const auto it_period = it_var_results->second.find(period);
                if (it_period == it_var_results->second.end()) continue;

                TGraphErrors* g =
                    make_kinematic_factor_graph(
                        it_period->second, period_color(period));
                g->SetMarkerSize(1.05);
                g->SetLineWidth(2);
                g->Draw("P SAME");
            }

            auto variable_label = std::make_unique<TLatex>();
            variable_label->SetNDC();
            variable_label->SetTextAlign(33);
            variable_label->SetTextFont(62);
            variable_label->SetTextSize(0.050);
            variable_label->DrawLatex(0.94, 0.93, v.title.c_str());
            owned.emplace_back(std::move(variable_label));
        }

        c.cd(8);
        draw_info_panel(sample_label, fit_version);

        c.cd(0);
        auto title = std::make_unique<TLatex>();
        title->SetNDC();
        title->SetTextAlign(22);
        title->SetTextFont(62);
        title->SetTextSize(0.028);
        title->DrawLatex(
            0.5, 0.990,
            (cfg.title + ": " + sample_label +
             " current-efficiency dependence on reconstructed kinematics").c_str());
        owned.emplace_back(std::move(title));

        c.Modified();
        c.Update();
        c.SaveAs((odir + "/" + output_name).c_str());
    };

    draw_one_canvas("c_kinematic_current_efficiency_data",
                    "current_efficiency_vs_kinematics.png",
                    "data",
                    data_results_by_var_period,
                    data_fits_by_var_period,
                    integrated_data,
                    false,
                    "Current-efficiency factor");

    draw_one_canvas("c_kinematic_current_efficiency_data_linear_fits",
                    "current_efficiency_vs_kinematics_linear_fits.png",
                    "data",
                    data_results_by_var_period,
                    data_fits_by_var_period,
                    integrated_data,
                    true,
                    "Current-efficiency factor");

    draw_one_canvas("c_kinematic_current_efficiency_mc",
                    "current_efficiency_vs_kinematics_mc.png",
                    "MC",
                    mc_results_by_var_period,
                    mc_fits_by_var_period,
                    integrated_mc,
                    false,
                    "Current-efficiency factor");

    draw_one_canvas("c_kinematic_current_efficiency_mc_linear_fits",
                    "current_efficiency_vs_kinematics_mc_linear_fits.png",
                    "MC",
                    mc_results_by_var_period,
                    mc_fits_by_var_period,
                    integrated_mc,
                    true,
                    "Current-efficiency factor");

    // Direct DATA/MC current-efficiency diagnostic.  This is the quantity
    // relevant to the cross-section acceptance correction and keeps the
    // presentation convention consistently DATA/MC.
    VarPeriodResults ratio_results_by_var_period;
    std::map<std::string, std::map<std::string, LinearFitSummary>> ratio_fits_by_var_period;
    std::map<std::string, double> integrated_ratio;

    for (const std::string& period : PERIOD_ORDER) {
        auto id = integrated_data.find(period);
        auto im = integrated_mc.find(period);
        if (id != integrated_data.end() && im != integrated_mc.end() &&
            std::isfinite(id->second) && id->second > 0.0 &&
            std::isfinite(im->second) && im->second > 0.0) {
            integrated_ratio[period] = id->second / im->second;
        }
    }

    std::ofstream ratio_csv(odir + "/current_efficiency_vs_kinematics_data_over_mc.csv");
    ratio_csv << "sample,period,variable,bin,x_low,x_high,x_center,data_over_mc,data_over_mc_stat,integrated_data_over_mc,n_data_current_points,n_mc_current_points,data_total_counts,mc_total_counts\n";

    std::ofstream ratio_fit_csv(odir + "/current_efficiency_vs_kinematics_data_over_mc_linear_fits.csv");
    ratio_fit_csv << "sample,period,variable,slope,intercept,slope_stat,intercept_stat,cov_slope_intercept,chi2,ndf,npoints\n";

    for (const KinematicVarConfig& v : vars) {
        for (const std::string& period : PERIOD_ORDER) {
            if (hide_sp19_inb_from_replacement_plots(hide_sp19_inb_from_all_period_plots, period)) {
                continue;
            }

            auto idv = data_results_by_var_period.find(v.key);
            auto imv = mc_results_by_var_period.find(v.key);
            if (idv == data_results_by_var_period.end() ||
                imv == mc_results_by_var_period.end()) continue;

            auto idp = idv->second.find(period);
            auto imp = imv->second.find(period);
            if (idp == idv->second.end() || imp == imv->second.end()) continue;

            std::vector<KinematicBinResult> rb =
                make_data_over_mc_kinematic_bins(idp->second, imp->second);
            ratio_results_by_var_period[v.key][period] = rb;

            const LinearFitSummary fit = fit_kinematic_factor_linear(rb);
            ratio_fits_by_var_period[v.key][period] = fit;
            ratio_fit_csv << "data_over_mc," << period << "," << v.key << ","
                          << fit.m << "," << fit.b << "," << fit.m_err << "," << fit.b_err << ","
                          << fit.cov_mb << "," << fit.chi2 << "," << fit.ndf << "," << fit.npoints << "\n";

            const size_t nb = std::min(idp->second.size(), imp->second.size());
            for (size_t ib = 0; ib < nb && ib < rb.size(); ++ib) {
                const KinematicBinResult& r = rb[ib];
                const KinematicBinResult& d = idp->second[ib];
                const KinematicBinResult& m = imp->second[ib];
                const double iratio = integrated_ratio.count(period)
                    ? integrated_ratio[period]
                    : std::numeric_limits<double>::quiet_NaN();
                ratio_csv << "data_over_mc," << period << "," << v.key << "," << ib << ","
                          << r.x_low << "," << r.x_high << "," << r.x_center << ","
                          << r.factor << "," << r.factor_err << "," << iratio << ","
                          << d.n_current_points << "," << m.n_current_points << ","
                          << d.total_counts << "," << m.total_counts << "\n";
            }
        }
    }

    draw_one_canvas("c_kinematic_current_efficiency_data_over_mc",
                    "current_efficiency_vs_kinematics_data_over_mc.png",
                    "DATA/MC",
                    ratio_results_by_var_period,
                    ratio_fits_by_var_period,
                    integrated_ratio,
                    false,
                    "DATA/MC current-efficiency ratio");

    draw_one_canvas("c_kinematic_current_efficiency_data_over_mc_linear_fits",
                    "current_efficiency_vs_kinematics_data_over_mc_linear_fits.png",
                    "DATA/MC",
                    ratio_results_by_var_period,
                    ratio_fits_by_var_period,
                    integrated_ratio,
                    true,
                    "DATA/MC current-efficiency ratio");

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

    int nth = std::max(1, std::min(7, max_workers));
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
        nth = std::max(1, std::min(7, max_workers));
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


// -----------------------------------------------------------------------------
// Photon-region current-dependence diagnostic (FT + six FD sectors)
// -----------------------------------------------------------------------------
//
// This study repeats the ep->epgamma current-dependence extraction in seven
// mutually exclusive photon regions.  These regional fits define the nominal
// event-level response model; the additional kinematic/theta studies remain
// diagnostic until explicitly promoted after validation.
//
//   <output_dir>/epg/sector_dependence_diagnostic/
//
// DATA/reconstructed MC region assignment uses detector2 and reconstructed
// p2_phi.  Generated MC has no reconstructed detector assignment, so the
// generated denominator is split geometrically using generated p2_theta and
// p2_phi: theta_gamma <= 5.5 deg is FT and larger angles are assigned to the
// standard six CLAS12 FD phi sectors.
// -----------------------------------------------------------------------------

static constexpr double PHOTON_REGION_FT_MAX_THETA_DEG = 5.5;

struct PhotonRegionSpec {
    std::string key;
    std::string label;
};

static const std::vector<PhotonRegionSpec> PHOTON_REGION_ORDER = {
    {"FT", "FT"},
    {"S1", "S1"},
    {"S2", "S2"},
    {"S3", "S3"},
    {"S4", "S4"},
    {"S5", "S5"},
    {"S6", "S6"}
};

static int photon_region_color(int i) {
    static const int colors[] = {
        kBlack,
        kBlue + 1,
        kOrange + 7,
        kGreen + 2,
        kRed + 1,
        kViolet + 1,
        kGray + 2
    };
    const int n = (int)(sizeof(colors) / sizeof(colors[0]));
    return colors[std::max(0, std::min(i, n - 1))];
}

static int photon_region_marker(int i) {
    static const int markers[] = {20, 21, 22, 23, 33, 34, 29};
    const int n = (int)(sizeof(markers) / sizeof(markers[0]));
    return markers[std::max(0, std::min(i, n - 1))];
}

static bool reconstructed_photon_region(const Branches& b,
                                        std::string& region_key) {
    if (!b.has_detector2) {
        return false;
    }

    if (b.detector2 == 0) {
        region_key = "FT";
        return true;
    }

    if (b.detector2 != 1 || !b.has_p2_phi) {
        return false;
    }

    const int sector = fd_sector_from_phi_rad(b.p2_phi);
    if (sector < 1 || sector > 6) {
        return false;
    }

    region_key = "S" + std::to_string(sector);
    return true;
}

static bool generated_photon_region(const Branches& b,
                                    std::string& region_key) {
    if (!(b.has_p2_theta && b.has_p2_phi)) {
        return false;
    }

    const double theta_deg = b.p2_theta * RAD2DEG;
    if (!std::isfinite(theta_deg) || !std::isfinite(b.p2_phi)) {
        return false;
    }

    if (theta_deg <= PHOTON_REGION_FT_MAX_THETA_DEG) {
        region_key = "FT";
        return true;
    }

    const int sector = fd_sector_from_phi_rad(b.p2_phi);
    if (sector < 1 || sector > 6) {
        return false;
    }

    region_key = "S" + std::to_string(sector);
    return true;
}

using PhotonRegionDataAggMap = std::map<std::string, DataAgg>;
using PhotonRegionMcCountMap = std::map<std::string, long long>;
using PhotonRegionResults = std::map<std::string, std::vector<PeriodResult>>;

static PhotonRegionDataAggMap process_data_tree_by_photon_region(
    const ChannelConfig& cfg,
    const std::string& key,
    TTree* tree,
    const std::unordered_map<int, ChargeEntry>& charge_map,
    const TopoCutMap& data_cuts,
    bool use_second_column_charge_for_all_unpolarized,
    bool use_columns_3_to_5_charge_sum_scaled_for_fa18_sp19_unpolarized,
    double columns_3_to_5_charge_sum_scale) {

    const PeriodTags tags = parse_period_from_key(key);
    PhotonRegionDataAggMap out;

    for (const PhotonRegionSpec& spec : PHOTON_REGION_ORDER) {
        out[spec.key].period = tags.display;
    }

    if (!tree) {
        return out;
    }

    Branches b;
    b.bind(tree);

    if (!b.has_runnum) {
        fatal("[current_dependence] FATAL: photon-region DATA diagnostic tree missing runnum.");
    }

    const Long64_t N = tree->GetEntries();
    std::unordered_map<int, ResolvedRunInfo> run_info_cache;
    run_info_cache.reserve(256);
    using ProgressClock = std::chrono::steady_clock;
    const auto progress_t0 = ProgressClock::now();
    const Long64_t progress_step = std::max<Long64_t>(250000, N / 10);
    Long64_t next_progress = progress_step;

    for (Long64_t i = 0; i < N; ++i) {
        tree->GetEntry(i);

        if (i + 1 >= next_progress || i + 1 == N) {
            const double sec = std::chrono::duration<double>(
                ProgressClock::now() - progress_t0).count();
            const double pct = N > 0 ? 100.0 * double(i + 1) / double(N) : 100.0;
            std::lock_guard<std::mutex> progress_lock(g_progress_output_mutex);
            std::cout << "[current_dependence] Photon-region DATA progress key="
                      << key << " " << (i + 1) << "/" << N
                      << " (" << std::fixed << std::setprecision(1) << pct
                      << "%) elapsed=" << sec << " s"
                      << std::defaultfloat << std::setprecision(6) << std::endl;
            next_progress += progress_step;
        }

        if (!passes_cone_cut(b)) {
            continue;
        }

        auto cache_it = run_info_cache.find(b.runnum);
        bool new_run = false;
        if (cache_it == run_info_cache.end()) {
            cache_it = run_info_cache.emplace(
                b.runnum,
                resolve_run_info_cached(
                    tags, b.runnum, charge_map,
                    use_second_column_charge_for_all_unpolarized,
                    use_columns_3_to_5_charge_sum_scaled_for_fa18_sp19_unpolarized,
                    columns_3_to_5_charge_sum_scale)).first;
            new_run = true;
        }

        const ResolvedRunInfo& run_info = cache_it->second;
        if (!run_info.valid) {
            continue;
        }

        if (new_run) {
            for (const PhotonRegionSpec& spec : PHOTON_REGION_ORDER) {
                DataAgg& agg = out[spec.key];
                agg.counts_by_run[b.runnum] += 0;
                agg.current_by_run[b.runnum] = run_info.current;
                agg.charge_by_run[b.runnum] = run_info.charge;
            }
        }

        if (!passes_global_dispatch(b, tags)) {
            continue;
        }

        if (!passes_sigma_dispatch(cfg, tags, data_cuts, b)) {
            continue;
        }

        std::string region_key;
        if (!reconstructed_photon_region(b, region_key)) {
            continue;
        }

        out[region_key].counts_by_run[b.runnum] += 1;
    }

    return out;
}

static PhotonRegionMcCountMap count_generated_tree_by_photon_region(const std::string& key, TTree* tree) {
    PhotonRegionMcCountMap counts;
    for (const PhotonRegionSpec& spec : PHOTON_REGION_ORDER) {
        counts[spec.key] = 0;
    }

    if (!tree) {
        return counts;
    }

    Branches b;
    b.bind(tree);

    if (!(b.has_p2_theta && b.has_p2_phi)) {
        fatal("[current_dependence] FATAL: generated MC photon-region diagnostic requires p2_theta and p2_phi.");
    }

    const Long64_t N = tree->GetEntries();
    using ProgressClock = std::chrono::steady_clock;
    const auto progress_t0 = ProgressClock::now();
    const Long64_t progress_step = std::max<Long64_t>(250000, N / 10);
    Long64_t next_progress = progress_step;

    for (Long64_t i = 0; i < N; ++i) {
        tree->GetEntry(i);

        if (i + 1 >= next_progress || i + 1 == N) {
            const double sec = std::chrono::duration<double>(
                ProgressClock::now() - progress_t0).count();
            const double pct = N > 0 ? 100.0 * double(i + 1) / double(N) : 100.0;
            std::lock_guard<std::mutex> progress_lock(g_progress_output_mutex);
            std::cout << "[current_dependence] Photon-region GEN progress key="
                      << key << " " << (i + 1) << "/" << N
                      << " (" << std::fixed << std::setprecision(1) << pct
                      << "%) elapsed=" << sec << " s"
                      << std::defaultfloat << std::setprecision(6) << std::endl;
            next_progress += progress_step;
        }

        std::string region_key;
        if (generated_photon_region(b, region_key)) {
            counts[region_key] += 1;
        }
    }

    return counts;
}

static PhotonRegionMcCountMap count_reconstructed_tree_by_photon_region(
    const ChannelConfig& cfg,
    const std::string& key,
    TTree* tree,
    const TopoCutMap& mc_cuts) {

    PhotonRegionMcCountMap counts;
    for (const PhotonRegionSpec& spec : PHOTON_REGION_ORDER) {
        counts[spec.key] = 0;
    }

    if (!tree) {
        return counts;
    }

    const PeriodTags tags = parse_period_from_key(key);
    Branches b;
    b.bind(tree);

    const Long64_t N = tree->GetEntries();
    using ProgressClock = std::chrono::steady_clock;
    const auto progress_t0 = ProgressClock::now();
    const Long64_t progress_step = std::max<Long64_t>(250000, N / 10);
    Long64_t next_progress = progress_step;

    for (Long64_t i = 0; i < N; ++i) {
        tree->GetEntry(i);

        if (i + 1 >= next_progress || i + 1 == N) {
            const double sec = std::chrono::duration<double>(
                ProgressClock::now() - progress_t0).count();
            const double pct = N > 0 ? 100.0 * double(i + 1) / double(N) : 100.0;
            std::lock_guard<std::mutex> progress_lock(g_progress_output_mutex);
            std::cout << "[current_dependence] Photon-region REC progress key="
                      << key << " " << (i + 1) << "/" << N
                      << " (" << std::fixed << std::setprecision(1) << pct
                      << "%) elapsed=" << sec << " s"
                      << std::defaultfloat << std::setprecision(6) << std::endl;
            next_progress += progress_step;
        }

        if (!passes_cone_cut(b)) {
            continue;
        }

        if (!passes_global_dispatch(b, tags)) {
            continue;
        }

        if (!passes_sigma_dispatch(cfg, tags, mc_cuts, b)) {
            continue;
        }

        std::string region_key;
        if (reconstructed_photon_region(b, region_key)) {
            counts[region_key] += 1;
        }
    }

    return counts;
}


using RegionThetaAggMap = std::map<std::string, std::map<std::string, std::vector<DataAgg>>>;
using PeriodRegionThetaAggMap = std::map<std::string, RegionThetaAggMap>;

static RegionThetaAggMap process_data_tree_region_theta(
    const ChannelConfig& cfg,
    const std::string& key,
    TTree* tree,
    const std::unordered_map<int, ChargeEntry>& charge_map,
    const TopoCutMap& data_cuts,
    const std::vector<KinematicVarConfig>& vars,
    bool use_second_column_charge_for_all_unpolarized,
    bool use_columns_3_to_5_charge_sum_scaled_for_fa18_sp19_unpolarized,
    double columns_3_to_5_charge_sum_scale) {

    RegionThetaAggMap out;
    const PeriodTags tags = parse_period_from_key(key);
    for (const PhotonRegionSpec& spec : PHOTON_REGION_ORDER) {
        for (const KinematicVarConfig& v : vars) {
            auto& bins = out[spec.key][v.key];
            bins.resize(v.edges.size() - 1);
            for (DataAgg& a : bins) a.period = tags.display;
        }
    }
    if (!tree || tags.display == "Fa18 Inb Supp") return out;

    Branches b;
    b.bind(tree);
    if (!b.has_runnum) fatal("[current_dependence] FATAL: region-theta DATA diagnostic tree missing runnum.");

    std::unordered_map<int, ResolvedRunInfo> run_info_cache;
    run_info_cache.reserve(256);
    const Long64_t N = tree->GetEntries();

    for (Long64_t i = 0; i < N; ++i) {
        tree->GetEntry(i);
        if (!passes_cone_cut(b)) continue;

        auto cache_it = run_info_cache.find(b.runnum);
        bool new_run = false;
        if (cache_it == run_info_cache.end()) {
            cache_it = run_info_cache.emplace(
                b.runnum,
                resolve_run_info_cached(tags, b.runnum, charge_map,
                    use_second_column_charge_for_all_unpolarized,
                    use_columns_3_to_5_charge_sum_scaled_for_fa18_sp19_unpolarized,
                    columns_3_to_5_charge_sum_scale)).first;
            new_run = true;
        }
        const ResolvedRunInfo& run_info = cache_it->second;
        if (!run_info.valid) continue;

        if (new_run) {
            for (const PhotonRegionSpec& spec : PHOTON_REGION_ORDER) {
                for (const KinematicVarConfig& v : vars) {
                    for (DataAgg& a : out[spec.key][v.key]) {
                        a.counts_by_run[b.runnum] += 0;
                        a.current_by_run[b.runnum] = run_info.current;
                        a.charge_by_run[b.runnum] = run_info.charge;
                    }
                }
            }
        }

        if (!passes_global_dispatch(b, tags)) continue;
        if (!passes_sigma_dispatch(cfg, tags, data_cuts, b)) continue;

        std::string region;
        if (!reconstructed_photon_region(b, region)) continue;

        for (const KinematicVarConfig& v : vars) {
            double value = std::numeric_limits<double>::quiet_NaN();
            if (!kinematic_value_for_config(b, v, value)) continue;
            const int ib = find_bin_index(v.edges, value);
            if (ib < 0) continue;
            DataAgg& a = out[region][v.key][ib];
            a.counts_by_run[b.runnum] += 1;
            a.kinematic_value_sum += value;
            a.kinematic_value_count += 1;
        }
    }
    return out;
}


using RegionThetaMcAggMap = std::map<std::string, KinematicMcAggMap>;
using PeriodRegionThetaMcAggMap = std::map<std::string, RegionThetaMcAggMap>;

static RegionThetaMcAggMap process_generated_tree_region_theta(
    const std::string& key,
    TTree* tree,
    const std::vector<KinematicVarConfig>& vars) {

    RegionThetaMcAggMap out;
    const PeriodTags tags = parse_period_from_key(key);
    const int current = parse_current_from_key(key);

    for (const PhotonRegionSpec& spec : PHOTON_REGION_ORDER) {
        for (const KinematicVarConfig& v : vars) {
            auto& bins = out[spec.key][v.key];
            bins.resize(v.edges.size() - 1);
            for (McKinematicBinAgg& a : bins) a.period = tags.display;
        }
    }

    if (!tree || tags.display == "Fa18 Inb Supp") return out;

    Branches b;
    b.bind(tree);
    const Long64_t N = tree->GetEntries();

    for (Long64_t i = 0; i < N; ++i) {
        tree->GetEntry(i);

        std::string region;
        if (!generated_photon_region(b, region)) continue;

        for (const KinematicVarConfig& v : vars) {
            double value = std::numeric_limits<double>::quiet_NaN();
            if (!kinematic_value_for_config(b, v, value)) continue;
            const int ib = find_bin_index(v.edges, value);
            if (ib < 0) continue;

            McKinematicBinAgg& bin = out[region][v.key][ib];
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

static RegionThetaMcAggMap process_reconstructed_tree_region_theta(
    const ChannelConfig& cfg,
    const std::string& key,
    TTree* tree,
    const TopoCutMap& mc_cuts,
    const std::vector<KinematicVarConfig>& vars) {

    RegionThetaMcAggMap out;
    const PeriodTags tags = parse_period_from_key(key);
    const int current = parse_current_from_key(key);

    for (const PhotonRegionSpec& spec : PHOTON_REGION_ORDER) {
        for (const KinematicVarConfig& v : vars) {
            auto& bins = out[spec.key][v.key];
            bins.resize(v.edges.size() - 1);
            for (McKinematicBinAgg& a : bins) a.period = tags.display;
        }
    }

    if (!tree || tags.display == "Fa18 Inb Supp") return out;

    Branches b;
    b.bind(tree);
    const Long64_t N = tree->GetEntries();

    for (Long64_t i = 0; i < N; ++i) {
        tree->GetEntry(i);
        if (!passes_cone_cut(b)) continue;
        if (!passes_global_dispatch(b, tags)) continue;
        if (!passes_sigma_dispatch(cfg, tags, mc_cuts, b)) continue;

        std::string region;
        if (!reconstructed_photon_region(b, region)) continue;

        for (const KinematicVarConfig& v : vars) {
            double value = std::numeric_limits<double>::quiet_NaN();
            if (!kinematic_value_for_config(b, v, value)) continue;
            const int ib = find_bin_index(v.edges, value);
            if (ib < 0) continue;

            McKinematicBinAgg& bin = out[region][v.key][ib];
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

static void merge_region_theta_mc(PeriodRegionThetaMcAggMap& merged,
                                  const std::string& period,
                                  const RegionThetaMcAggMap& one,
                                  const std::vector<KinematicVarConfig>& vars) {
    for (const PhotonRegionSpec& spec : PHOTON_REGION_ORDER) {
        auto ir = one.find(spec.key);
        if (ir == one.end()) continue;
        for (const KinematicVarConfig& v : vars) {
            auto iv = ir->second.find(v.key);
            if (iv == ir->second.end()) continue;
            auto& dst = merged[period][spec.key][v.key];
            const auto& src = iv->second;
            const size_t n = std::min(dst.size(), src.size());
            for (size_t ib = 0; ib < n; ++ib) merge_mc_kinematic_bin(dst[ib], src[ib]);
        }
    }
}

struct RelativeSlopePoint {
    int region_index = -1;
    double x = std::numeric_limits<double>::quiet_NaN();
    double slope = std::numeric_limits<double>::quiet_NaN();
    double slope_err = std::numeric_limits<double>::quiet_NaN();
    int n_current_points = 0;
    double total_counts = 0.0;
};

struct PooledRegionThetaSlopeFit {
    bool valid = false;
    std::array<double, 7> region_intercepts{};
    std::array<double, 7> region_intercept_err{};
    std::array<double, 7> region_centers{};
    double theta_gradient = std::numeric_limits<double>::quiet_NaN();
    double theta_gradient_err = std::numeric_limits<double>::quiet_NaN();
    double chi2 = std::numeric_limits<double>::quiet_NaN();
    int ndf = 0;
    int npoints = 0;
};

using PeriodPooledAngleFitMap = std::map<std::string, PooledRegionThetaSlopeFit>;

static bool invert_small_matrix(std::vector<std::vector<double>> a,
                                std::vector<std::vector<double>>& inv) {
    const int n = (int)a.size();
    inv.assign(n, std::vector<double>(n, 0.0));
    for (int i = 0; i < n; ++i) inv[i][i] = 1.0;

    for (int col = 0; col < n; ++col) {
        int pivot = col;
        for (int r = col + 1; r < n; ++r)
            if (std::fabs(a[r][col]) > std::fabs(a[pivot][col])) pivot = r;
        if (!(std::fabs(a[pivot][col]) > 1e-18) || !std::isfinite(a[pivot][col])) return false;
        if (pivot != col) {
            std::swap(a[pivot], a[col]);
            std::swap(inv[pivot], inv[col]);
        }

        const double p = a[col][col];
        for (int j = 0; j < n; ++j) {
            a[col][j] /= p;
            inv[col][j] /= p;
        }
        for (int r = 0; r < n; ++r) {
            if (r == col) continue;
            const double f = a[r][col];
            for (int j = 0; j < n; ++j) {
                a[r][j] -= f * a[col][j];
                inv[r][j] -= f * inv[col][j];
            }
        }
    }
    return true;
}

static PooledRegionThetaSlopeFit fit_pooled_region_theta_relative_slopes(
    const std::vector<RelativeSlopePoint>& pts) {

    PooledRegionThetaSlopeFit out;
    constexpr int NP = 8; // seven region intercepts + one common theta gradient

    std::array<double, 7> sw{};
    std::array<double, 7> swx{};
    for (const RelativeSlopePoint& p : pts) {
        if (p.region_index < 0 || p.region_index >= 7 ||
            !std::isfinite(p.x) || !std::isfinite(p.slope) ||
            !std::isfinite(p.slope_err) || !(p.slope_err > 0.0)) continue;
        const double w = 1.0 / (p.slope_err * p.slope_err);
        sw[p.region_index] += w;
        swx[p.region_index] += w * p.x;
    }
    for (int r = 0; r < 7; ++r) {
        out.region_centers[r] = (sw[r] > 0.0) ? swx[r] / sw[r] : 0.0;
    }

    std::vector<std::vector<double>> normal(NP, std::vector<double>(NP, 0.0));
    std::vector<double> rhs(NP, 0.0);
    std::vector<RelativeSlopePoint> used;

    for (const RelativeSlopePoint& p : pts) {
        if (p.region_index < 0 || p.region_index >= 7 ||
            !std::isfinite(p.x) || !std::isfinite(p.slope) ||
            !std::isfinite(p.slope_err) || !(p.slope_err > 0.0)) continue;

        const double w = 1.0 / (p.slope_err * p.slope_err);
        std::array<double, NP> row{};
        row[p.region_index] = 1.0;
        row[7] = p.x - out.region_centers[p.region_index];

        for (int i = 0; i < NP; ++i) {
            rhs[i] += w * row[i] * p.slope;
            for (int j = 0; j < NP; ++j) normal[i][j] += w * row[i] * row[j];
        }
        used.push_back(p);
    }

    out.npoints = (int)used.size();
    if (out.npoints <= NP) return out;
    for (int r = 0; r < 7; ++r) if (!(sw[r] > 0.0)) return out;

    std::vector<std::vector<double>> cov;
    if (!invert_small_matrix(normal, cov)) return out;

    std::vector<double> beta(NP, 0.0);
    for (int i = 0; i < NP; ++i)
        for (int j = 0; j < NP; ++j)
            beta[i] += cov[i][j] * rhs[j];

    for (int r = 0; r < 7; ++r) {
        out.region_intercepts[r] = beta[r];
        out.region_intercept_err[r] = std::sqrt(std::max(0.0, cov[r][r]));
    }
    out.theta_gradient = beta[7];
    out.theta_gradient_err = std::sqrt(std::max(0.0, cov[7][7]));

    double chi2 = 0.0;
    for (const RelativeSlopePoint& p : used) {
        const double pred =
            out.region_intercepts[p.region_index] +
            out.theta_gradient * (p.x - out.region_centers[p.region_index]);
        const double pull = (p.slope - pred) / p.slope_err;
        chi2 += pull * pull;
    }
    out.chi2 = chi2;
    out.ndf = out.npoints - NP;
    out.valid = std::isfinite(out.theta_gradient) &&
                std::isfinite(out.theta_gradient_err) &&
                out.theta_gradient_err >= 0.0;
    return out;
}


struct PooledRegionQuadraticFit {
    bool valid = false;
    std::array<double, 7> region_intercepts{};
    std::array<double, 7> region_centers{};
    std::array<double, 7> region_second_moments{};
    double linear = std::numeric_limits<double>::quiet_NaN();
    double linear_err = std::numeric_limits<double>::quiet_NaN();
    double quadratic = std::numeric_limits<double>::quiet_NaN();
    double quadratic_err = std::numeric_limits<double>::quiet_NaN();
    double chi2 = std::numeric_limits<double>::quiet_NaN();
    int ndf = 0;
    int npoints = 0;
    double aic = std::numeric_limits<double>::quiet_NaN();
    double bic = std::numeric_limits<double>::quiet_NaN();
};

static PooledRegionQuadraticFit fit_pooled_region_quadratic_relative_slopes(
    const std::vector<RelativeSlopePoint>& pts) {

    PooledRegionQuadraticFit out;
    constexpr int NP = 9; // seven regional baselines + common linear + common quadratic

    std::array<double, 7> sw{};
    std::array<double, 7> swx{};
    std::vector<RelativeSlopePoint> used;
    used.reserve(pts.size());

    for (const RelativeSlopePoint& p : pts) {
        if (p.region_index < 0 || p.region_index >= 7 ||
            !std::isfinite(p.x) || !std::isfinite(p.slope) ||
            !std::isfinite(p.slope_err) || !(p.slope_err > 0.0)) continue;
        const double w = 1.0 / (p.slope_err * p.slope_err);
        sw[p.region_index] += w;
        swx[p.region_index] += w * p.x;
        used.push_back(p);
    }

    out.npoints = (int)used.size();
    if (out.npoints <= NP) return out;
    for (int r = 0; r < 7; ++r) {
        if (!(sw[r] > 0.0)) return out;
        out.region_centers[r] = swx[r] / sw[r];
    }

    // Center both the linear and quadratic basis within each photon region.
    // This keeps the regional intercept equal to the weighted regional mean
    // response while the added terms describe only residual kinematic shape.
    std::array<double, 7> swz2{};
    for (const RelativeSlopePoint& p : used) {
        const double w = 1.0 / (p.slope_err * p.slope_err);
        const double z = p.x - out.region_centers[p.region_index];
        swz2[p.region_index] += w * z * z;
    }
    for (int r = 0; r < 7; ++r)
        out.region_second_moments[r] = swz2[r] / sw[r];

    std::vector<std::vector<double>> normal(NP, std::vector<double>(NP, 0.0));
    std::vector<double> rhs(NP, 0.0);

    for (const RelativeSlopePoint& p : used) {
        const double w = 1.0 / (p.slope_err * p.slope_err);
        const int r = p.region_index;
        const double z = p.x - out.region_centers[r];
        const double q = z * z - out.region_second_moments[r];

        std::array<double, NP> row{};
        row[r] = 1.0;
        row[7] = z;
        row[8] = q;

        for (int i = 0; i < NP; ++i) {
            rhs[i] += w * row[i] * p.slope;
            for (int j = 0; j < NP; ++j)
                normal[i][j] += w * row[i] * row[j];
        }
    }

    std::vector<std::vector<double>> cov;
    if (!invert_small_matrix(normal, cov)) return out;

    std::vector<double> beta(NP, 0.0);
    for (int i = 0; i < NP; ++i)
        for (int j = 0; j < NP; ++j)
            beta[i] += cov[i][j] * rhs[j];

    for (int r = 0; r < 7; ++r)
        out.region_intercepts[r] = beta[r];
    out.linear = beta[7];
    out.linear_err = std::sqrt(std::max(0.0, cov[7][7]));
    out.quadratic = beta[8];
    out.quadratic_err = std::sqrt(std::max(0.0, cov[8][8]));

    out.chi2 = 0.0;
    for (const RelativeSlopePoint& p : used) {
        const int r = p.region_index;
        const double z = p.x - out.region_centers[r];
        const double q = z * z - out.region_second_moments[r];
        const double pred =
            out.region_intercepts[r] + out.linear * z + out.quadratic * q;
        const double pull = (p.slope - pred) / p.slope_err;
        out.chi2 += pull * pull;
    }

    out.ndf = out.npoints - NP;
    out.aic = out.chi2 + 2.0 * NP;
    out.bic = out.chi2 + double(NP) * std::log(double(out.npoints));
    out.valid = std::isfinite(out.chi2) &&
                std::isfinite(out.linear) &&
                std::isfinite(out.quadratic);
    return out;
}



struct RelativeSlopeModelComparison {
    bool valid = false;
    double chi2_m0 = std::numeric_limits<double>::quiet_NaN();
    int ndf_m0 = 0;
    double aic_m0 = std::numeric_limits<double>::quiet_NaN();
    double bic_m0 = std::numeric_limits<double>::quiet_NaN();

    double chi2_m1 = std::numeric_limits<double>::quiet_NaN();
    int ndf_m1 = 0;
    double aic_m1 = std::numeric_limits<double>::quiet_NaN();
    double bic_m1 = std::numeric_limits<double>::quiet_NaN();

    double chi2_m2 = std::numeric_limits<double>::quiet_NaN();
    int ndf_m2 = 0;
    double aic_m2 = std::numeric_limits<double>::quiet_NaN();
    double bic_m2 = std::numeric_limits<double>::quiet_NaN();

    double gradient = std::numeric_limits<double>::quiet_NaN();
    double gradient_err = std::numeric_limits<double>::quiet_NaN();
    double gradient_significance = std::numeric_limits<double>::quiet_NaN();
    int npoints = 0;
};

static RelativeSlopeModelComparison compare_relative_slope_models(
    const std::vector<RelativeSlopePoint>& pts) {

    RelativeSlopeModelComparison out;
    std::vector<RelativeSlopePoint> used;
    used.reserve(pts.size());
    for (const RelativeSlopePoint& p : pts) {
        if (p.region_index < 0 || p.region_index >= 7 ||
            !std::isfinite(p.x) || !std::isfinite(p.slope) ||
            !std::isfinite(p.slope_err) || !(p.slope_err > 0.0)) continue;
        used.push_back(p);
    }
    out.npoints = (int)used.size();
    if (out.npoints <= 8) return out;

    auto info_criteria = [&](double chi2, int k, double& aic, double& bic) {
        aic = chi2 + 2.0 * double(k);
        bic = chi2 + double(k) * std::log(double(out.npoints));
    };

    // M0: one common relative current slope for the entire period.
    double sw = 0.0;
    double swy = 0.0;
    for (const RelativeSlopePoint& p : used) {
        const double w = 1.0 / (p.slope_err * p.slope_err);
        sw += w;
        swy += w * p.slope;
    }
    if (!(sw > 0.0)) return out;
    const double common = swy / sw;
    out.chi2_m0 = 0.0;
    for (const RelativeSlopePoint& p : used) {
        const double pull = (p.slope - common) / p.slope_err;
        out.chi2_m0 += pull * pull;
    }
    out.ndf_m0 = out.npoints - 1;
    info_criteria(out.chi2_m0, 1, out.aic_m0, out.bic_m0);

    // M1: one constant relative current slope for each photon region.
    std::array<double, 7> region_sw{};
    std::array<double, 7> region_swy{};
    for (const RelativeSlopePoint& p : used) {
        const double w = 1.0 / (p.slope_err * p.slope_err);
        region_sw[p.region_index] += w;
        region_swy[p.region_index] += w * p.slope;
    }
    int populated_regions = 0;
    std::array<double, 7> region_mean{};
    for (int ir = 0; ir < 7; ++ir) {
        if (region_sw[ir] > 0.0) {
            region_mean[ir] = region_swy[ir] / region_sw[ir];
            ++populated_regions;
        }
    }
    if (populated_regions != 7) return out;

    out.chi2_m1 = 0.0;
    for (const RelativeSlopePoint& p : used) {
        const double pull = (p.slope - region_mean[p.region_index]) / p.slope_err;
        out.chi2_m1 += pull * pull;
    }
    out.ndf_m1 = out.npoints - 7;
    info_criteria(out.chi2_m1, 7, out.aic_m1, out.bic_m1);

    // M2: seven regional baselines plus one common linear dependence on x.
    const PooledRegionThetaSlopeFit pooled =
        fit_pooled_region_theta_relative_slopes(used);
    if (!pooled.valid) return out;

    out.chi2_m2 = pooled.chi2;
    out.ndf_m2 = pooled.ndf;
    info_criteria(out.chi2_m2, 8, out.aic_m2, out.bic_m2);
    out.gradient = pooled.theta_gradient;
    out.gradient_err = pooled.theta_gradient_err;
    if (std::isfinite(out.gradient) && std::isfinite(out.gradient_err) &&
        out.gradient_err > 0.0) {
        out.gradient_significance = out.gradient / out.gradient_err;
    }
    out.valid = true;
    return out;
}

static std::vector<RelativeSlopePoint> relative_slope_points_from_data_theta_bins(
    const std::map<std::string, std::vector<DataAgg>>& by_region,
    const KinematicVarConfig& theta_var) {

    std::vector<RelativeSlopePoint> out;
    for (int ir = 0; ir < 7; ++ir) {
        auto it = by_region.find(PHOTON_REGION_ORDER[ir].key);
        if (it == by_region.end()) continue;
        const auto& bins = it->second;
        for (size_t ib = 0; ib < bins.size() && ib + 1 < theta_var.edges.size(); ++ib) {
            const auto pts = data_points_from_agg(bins[ib]);
            if (pts.size() < 2) continue;
            const FitResult fit = fit_points(pts);
            const double s = fit_percent_slope(fit) / 100.0;
            const double se = fit_percent_slope_err(fit) / 100.0;
            if (!std::isfinite(s) || !std::isfinite(se) || !(se > 0.0)) continue;

            RelativeSlopePoint p;
            p.region_index = ir;
            p.x = mean_or_bin_center(theta_var.edges[ib], theta_var.edges[ib+1],
                                     bins[ib].kinematic_value_sum,
                                     bins[ib].kinematic_value_count);
            p.slope = s;
            p.slope_err = se;
            p.n_current_points = (int)pts.size();
            p.total_counts = (double)bins[ib].kinematic_value_count;
            out.push_back(p);
        }
    }
    return out;
}

static std::vector<RelativeSlopePoint> relative_slope_points_from_mc_theta_bins(
    const std::map<std::string, std::vector<McKinematicBinAgg>>& by_region,
    const KinematicVarConfig& theta_var) {

    std::vector<RelativeSlopePoint> out;
    for (int ir = 0; ir < 7; ++ir) {
        auto it = by_region.find(PHOTON_REGION_ORDER[ir].key);
        if (it == by_region.end()) continue;
        const auto& bins = it->second;
        for (size_t ib = 0; ib < bins.size() && ib + 1 < theta_var.edges.size(); ++ib) {
            std::vector<McAgg> aggs;
            for (const auto& kv : bins[ib].by_current) aggs.push_back(kv.second);
            const auto pts = mc_points_from_aggs(aggs, bins[ib].period);
            if (pts.size() < 2) continue;
            const FitResult fit = fit_points(pts);
            const double s = fit_percent_slope(fit) / 100.0;
            const double se = fit_percent_slope_err(fit) / 100.0;
            if (!std::isfinite(s) || !std::isfinite(se) || !(se > 0.0)) continue;

            RelativeSlopePoint p;
            p.region_index = ir;
            p.x = bins[ib].generated_kinematic_value_count > 0
                ? bins[ib].generated_kinematic_value_sum / double(bins[ib].generated_kinematic_value_count)
                : mean_or_bin_center(theta_var.edges[ib], theta_var.edges[ib+1],
                                     bins[ib].reconstructed_kinematic_value_sum,
                                     bins[ib].reconstructed_kinematic_value_count);
            p.slope = s;
            p.slope_err = se;
            p.n_current_points = (int)pts.size();
            double total_rec = 0.0;
            for (const auto& kv : bins[ib].by_current)
                total_rec += (double)kv.second.n_rec;
            p.total_counts = total_rec;
            out.push_back(p);
        }
    }
    return out;
}

static PeriodPooledAngleFitMap run_region_theta_data_diagnostic(
    const ChannelConfig& cfg,
    const std::map<std::string, TTree*>& data_trees,
    const std::map<std::string, TTree*>& gen_trees,
    const std::map<std::string, TTree*>& rec_trees,
    const std::unordered_map<int, ChargeEntry>& charge_map,
    const TopoCutMap& data_cuts,
    const TopoCutMap& mc_cuts,
    const std::string& output_dir,
    int max_workers,
    bool use_second_column_charge_for_all_unpolarized,
    bool use_columns_3_to_5_charge_sum_scaled_for_fa18_sp19_unpolarized,
    double columns_3_to_5_charge_sum_scale,
    bool hide_sp19,
    bool process_mc) {

    (void)hide_sp19;
    // Use the same complete kinematic set as the inclusive current-dependence
    // diagnostic.  The region-conditioned scan is used for model selection:
    // it tests whether a residual dependence remains after the photon-region
    // response is accounted for, without assuming in advance that theta_e or
    // theta_gamma is the relevant coordinate.
    const std::vector<KinematicVarConfig> vars = kinematic_current_var_configs();
    PeriodPooledAngleFitMap production_e_theta_models;

    PeriodRegionThetaAggMap merged;
    for (const std::string& period : PERIOD_ORDER) {
        for (const PhotonRegionSpec& spec : PHOTON_REGION_ORDER) {
            for (const KinematicVarConfig& v : vars) {
                auto& bins = merged[period][spec.key][v.key];
                bins.resize(v.edges.size() - 1);
                for (DataAgg& a : bins) a.period = period;
            }
        }
    }

    std::vector<std::pair<std::string, TTree*>> items;
    for (const auto& kv : data_trees) {
        const PeriodTags tags = parse_period_from_key(kv.first);
        if (tags.display != "Fa18 Inb Supp") items.push_back(kv);
    }
    std::mutex merge_mutex;
    std::mutex progress_mutex;
    std::atomic<int> done{0};
    int nth = std::max(1, std::min(7, max_workers));
    nth = std::min(nth, std::max(1, (int)items.size()));
    using Clock = std::chrono::steady_clock;
    const auto phase_t0 = Clock::now();

    std::cout << "[current_dependence] Region-theta DATA diagnostic for "
              << cfg.csv_channel << ": " << items.size()
              << " tree(s), " << nth << " worker(s), variables=Q2,xB,t,phi,e_theta,p_theta,g_theta."
              << std::endl;

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(nth)
#endif
    for (int i = 0; i < (int)items.size(); ++i) {
        const PeriodTags tags = parse_period_from_key(items[i].first);
        const auto tree_t0 = Clock::now();
        const Long64_t nentries = items[i].second ? items[i].second->GetEntries() : 0;
        {
            std::lock_guard<std::mutex> lock(progress_mutex);
            std::cout << "[current_dependence] Region-theta DATA start "
                      << (i + 1) << "/" << items.size()
                      << " key=" << items[i].first
                      << " entries=" << nentries << std::endl;
        }

        RegionThetaAggMap local = process_data_tree_region_theta(
            cfg, items[i].first, items[i].second, charge_map, data_cuts, vars,
            use_second_column_charge_for_all_unpolarized,
            use_columns_3_to_5_charge_sum_scaled_for_fa18_sp19_unpolarized,
            columns_3_to_5_charge_sum_scale);

        {
            std::lock_guard<std::mutex> lock(merge_mutex);
            for (const PhotonRegionSpec& spec : PHOTON_REGION_ORDER) {
                for (const KinematicVarConfig& v : vars) {
                    auto& dst = merged[tags.display][spec.key][v.key];
                    const auto& src = local[spec.key][v.key];
                    for (size_t ib = 0; ib < dst.size() && ib < src.size(); ++ib) {
                        merge_data_agg(dst[ib], src[ib]);
                    }
                }
            }
        }

        const int finished = ++done;
        {
            std::lock_guard<std::mutex> lock(progress_mutex);
            const double sec = std::chrono::duration<double>(Clock::now() - tree_t0).count();
            std::cout << "[current_dependence] Region-theta DATA done "
                      << finished << "/" << items.size()
                      << " key=" << items[i].first
                      << " elapsed=" << std::fixed << std::setprecision(1)
                      << sec << " s" << std::defaultfloat << std::setprecision(6) << std::endl;
        }
    }

    std::cout << "[current_dependence] Region-theta DATA scan complete for "
              << cfg.csv_channel << " in "
              << std::fixed << std::setprecision(1)
              << std::chrono::duration<double>(Clock::now() - phase_t0).count()
              << " s." << std::defaultfloat << std::setprecision(6) << std::endl;

    PeriodRegionThetaMcAggMap merged_mc;
    if (process_mc) {
        for (const std::string& period : PERIOD_ORDER) {
            for (const PhotonRegionSpec& spec : PHOTON_REGION_ORDER) {
                for (const KinematicVarConfig& v : vars) {
                    auto& bins = merged_mc[period][spec.key][v.key];
                    bins.resize(v.edges.size() - 1);
                    for (McKinematicBinAgg& a : bins) a.period = period;
                }
            }
        }

        std::mutex mc_merge_mutex;

        std::vector<std::pair<std::string, TTree*>> gen_items;
        for (const auto& kv : gen_trees) {
            const PeriodTags tags = parse_period_from_key(kv.first);
            if (tags.display != "Fa18 Inb Supp") gen_items.push_back(kv);
        }
        int mc_nth = std::max(1, std::min(7, max_workers));
        mc_nth = std::min(mc_nth, std::max(1, (int)gen_items.size()));
        std::cout << "[current_dependence] Region-theta generated-MC diagnostic for "
                  << cfg.csv_channel << ": " << gen_items.size()
                  << " tree(s), " << mc_nth << " worker(s)." << std::endl;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(mc_nth)
#endif
        for (int i = 0; i < (int)gen_items.size(); ++i) {
            const PeriodTags tags = parse_period_from_key(gen_items[i].first);
            RegionThetaMcAggMap one =
                process_generated_tree_region_theta(gen_items[i].first, gen_items[i].second, vars);
            std::lock_guard<std::mutex> lock(mc_merge_mutex);
            merge_region_theta_mc(merged_mc, tags.display, one, vars);
        }

        std::vector<std::pair<std::string, TTree*>> rec_items;
        for (const auto& kv : rec_trees) {
            const PeriodTags tags = parse_period_from_key(kv.first);
            if (tags.display != "Fa18 Inb Supp") rec_items.push_back(kv);
        }
        mc_nth = std::max(1, std::min(7, max_workers));
        mc_nth = std::min(mc_nth, std::max(1, (int)rec_items.size()));
        std::cout << "[current_dependence] Region-theta reconstructed-MC diagnostic for "
                  << cfg.csv_channel << ": " << rec_items.size()
                  << " tree(s), " << mc_nth << " worker(s)." << std::endl;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(mc_nth)
#endif
        for (int i = 0; i < (int)rec_items.size(); ++i) {
            const PeriodTags tags = parse_period_from_key(rec_items[i].first);
            RegionThetaMcAggMap one =
                process_reconstructed_tree_region_theta(
                    cfg, rec_items[i].first, rec_items[i].second, mc_cuts, vars);
            std::lock_guard<std::mutex> lock(mc_merge_mutex);
            merge_region_theta_mc(merged_mc, tags.display, one, vars);
        }
    }

    const std::string odir = output_dir + "/" + cfg.output_token + "/sector_dependence_diagnostic/region_theta";
    mkdir_p(odir);

    std::ofstream data_csv(odir + "/data_current_efficiency_by_region_and_theta.csv");
    data_csv << "period,region,variable,bin,x_low,x_high,x_center,factor,factor_stat,n_current_points,total_counts\n";

    std::ofstream mc_csv;
    std::ofstream ratio_csv;
    if (process_mc) {
        mc_csv.open(odir + "/mc_current_efficiency_by_region_and_theta.csv");
        mc_csv << "period,region,variable,bin,x_low,x_high,x_center,factor,factor_stat,n_current_points,total_counts\n";
        ratio_csv.open(odir + "/data_over_mc_current_efficiency_by_region_and_theta.csv");
        ratio_csv << "period,region,variable,bin,x_low,x_high,x_center,data_factor,data_factor_stat,mc_factor,mc_factor_stat,data_over_mc,data_over_mc_stat,data_n_current_points,mc_n_current_points\n";
    }

    const std::vector<std::string> periods = {"Sp18 Inb", "Sp18 Out", "Fa18 Inb", "Fa18 Out"};
    const std::string note_dir = output_dir + "/analysis_note";
    mkdir_p(note_dir);

    auto draw_region_theta_canvas =
        [&](const KinematicVarConfig& v,
            const std::string& sample_token,
            const std::string& y_title,
            const std::map<std::string, std::map<std::string, std::vector<KinematicBinResult>>>& results,
            double fallback_ymin,
            double fallback_ymax) {

        // One fixed production scale is used for all run periods.  Large
        // low-statistics error bars are intentionally allowed to clip rather
        // than stretching the full canvas.
        const auto yrange = production_efficiency_plot_range(
            sample_token == "data_over_mc" ? "data_over_mc" : "factor");
        const double ymin = yrange.first;
        const double ymax = yrange.second;

        TCanvas c(("c_region_theta_" + cfg.output_token + "_" + sample_token + "_" + v.key).c_str(),
                  "", 1800, 1200);
        c.Divide(2, 2, 0.002, 0.002);
        std::vector<std::unique_ptr<TObject>> owned;

        for (size_t ip = 0; ip < periods.size(); ++ip) {
            c.cd((int)ip + 1);
            configure_production_pad((TPad*)gPad, 0.145, 0.035, 0.145, 0.075);

            auto frame = std::make_unique<TH1D>(
                ("h_rt_" + cfg.output_token + sample_token + v.key + std::to_string(ip)).c_str(),
                (";" + v.x_label + ";" + y_title).c_str(),
                100, v.edges.front(), v.edges.back());
            frame->SetDirectory(nullptr);
            frame->SetMinimum(ymin);
            frame->SetMaximum(ymax);
            style_production_frame(frame.get(), 0.050, 0.050, 0.041);
            frame->Draw();
            owned.emplace_back(std::move(frame));

            auto unity = std::make_unique<TLine>(v.edges.front(), 1.0, v.edges.back(), 1.0);
            unity->SetLineStyle(2);
            unity->SetLineWidth(2);
            unity->SetLineColor(kGray + 1);
            unity->Draw("SAME");
            owned.emplace_back(std::move(unity));

            auto iper = results.find(periods[ip]);
            std::vector<TGraphErrors*> legend_graphs;
            if (iper != results.end()) {
                for (int ir = 0; ir < 7; ++ir) {
                    auto ireg = iper->second.find(PHOTON_REGION_ORDER[ir].key);
                    if (ireg == iper->second.end()) continue;
                    std::unique_ptr<TGraphErrors> graph(
                        make_region_theta_factor_graph(ireg->second, photon_region_color(ir)));
                    graph->SetMarkerStyle(photon_region_marker(ir));
                    graph->SetMarkerSize(1.15);
                    graph->SetLineWidth(2);
                    graph->Draw("P SAME");
                    if (ip == 0) legend_graphs.push_back(graph.get());
                    owned.emplace_back(std::move(graph));
                }
            }

            // One shared mapping legend is sufficient because marker/color
            // conventions are identical in all four panels.
            if (ip == 0) {
                auto leg = std::make_unique<TLegend>(0.17, 0.73, 0.78, 0.89);
                leg->SetNColumns(4);
                leg->SetBorderSize(0);
                leg->SetFillStyle(0);
                leg->SetTextFont(42);
                leg->SetTextSize(0.032);
                for (size_t ir = 0; ir < legend_graphs.size() && ir < PHOTON_REGION_ORDER.size(); ++ir)
                    leg->AddEntry(legend_graphs[ir], PHOTON_REGION_ORDER[ir].label.c_str(), "pe");
                leg->Draw();
                owned.emplace_back(std::move(leg));
            }

            owned.emplace_back(make_period_label(periods[ip], 0.94, 0.86, 33, 0.052));
        }

        c.cd(0);
        auto title = std::make_unique<TLatex>();
        title->SetNDC();
        title->SetTextAlign(22);
        title->SetTextFont(62);
        title->SetTextSize(0.030);
        std::string sample_title = sample_token;
        if (sample_token == "data") sample_title = "DATA";
        else if (sample_token == "mc") sample_title = "MC";
        else if (sample_token == "data_over_mc") sample_title = "DATA/MC";
        title->DrawLatex(
            0.5, 0.972,
            (cfg.title + ": " + sample_title + " current efficiency vs. " +
             v.title + " within photon regions").c_str());
        owned.emplace_back(std::move(title));

        c.Modified();
        c.Update();
        const std::string path =
            odir + "/" + sample_token + "_current_efficiency_by_region_vs_" + v.key + ".png";
        c.SaveAs(path.c_str());

        const std::string note_path =
            note_dir + "/" + cfg.output_token + "_regional_" + sample_token +
            "_current_efficiency_vs_" + v.key + ".png";
        gSystem->CopyFile(path.c_str(), note_path.c_str(), true);
    };

    std::ofstream model_csv(odir + "/regional_kinematic_model_comparison.csv");
    model_csv
        << "sample,period,variable,npoints,"
        << "m0_chi2,m0_ndf,m0_aic,m0_bic,"
        << "m1_chi2,m1_ndf,m1_aic,m1_bic,"
        << "m2_chi2,m2_ndf,m2_aic,m2_bic,"
        << "gradient,gradient_stat,gradient_significance,"
        << "delta_chi2_m1_to_m2,delta_aic_m1_to_m2,delta_bic_m1_to_m2\n";

    struct CandidateModelRow {
        std::string sample;
        std::string period;
        std::string variable;
        RelativeSlopeModelComparison cmp;
    };
    std::vector<CandidateModelRow> model_rows;

    using PeriodRegionResultMap =
        std::map<std::string, std::map<std::string, std::vector<KinematicBinResult>>>;
    std::map<std::string, PeriodRegionResultMap> stored_data_results;
    std::map<std::string, PeriodRegionResultMap> stored_mc_results;
    std::map<std::string, PeriodRegionResultMap> stored_ratio_results;

    for (const KinematicVarConfig& v : vars) {
        std::map<std::string, std::map<std::string, std::vector<KinematicBinResult>>> data_results;
        std::map<std::string, std::map<std::string, std::vector<KinematicBinResult>>> mc_results;
        std::map<std::string, std::map<std::string, std::vector<KinematicBinResult>>> ratio_results;

        for (const std::string& period : periods) {
            for (int ir = 0; ir < 7; ++ir) {
                const std::string& rkey = PHOTON_REGION_ORDER[ir].key;
                const auto& dbins = merged[period][rkey][v.key];
                std::vector<KinematicBinResult> dres;
                for (size_t ib = 0; ib < dbins.size(); ++ib) {
                    KinematicBinResult br =
                        current_factor_for_kinematic_bin(dbins[ib], v.edges[ib], v.edges[ib+1]);
                    dres.push_back(br);
                    data_csv << period << "," << rkey << "," << v.key << "," << ib << ","
                             << std::setprecision(12)
                             << br.x_low << "," << br.x_high << "," << br.x_center << ","
                             << br.factor << "," << br.factor_err << ","
                             << br.n_current_points << "," << br.total_counts << "\n";
                }
                data_results[period][rkey] = dres;

                if (process_mc) {
                    const auto& mbins = merged_mc[period][rkey][v.key];
                    std::vector<KinematicBinResult> mres;
                    for (size_t ib = 0; ib < mbins.size(); ++ib) {
                        KinematicBinResult br =
                            current_factor_for_kinematic_mc_bin(mbins[ib], v.edges[ib], v.edges[ib+1]);
                        mres.push_back(br);
                        mc_csv << period << "," << rkey << "," << v.key << "," << ib << ","
                               << std::setprecision(12)
                               << br.x_low << "," << br.x_high << "," << br.x_center << ","
                               << br.factor << "," << br.factor_err << ","
                               << br.n_current_points << "," << br.total_counts << "\n";
                    }
                    mc_results[period][rkey] = mres;

                    const auto rres = make_data_over_mc_kinematic_bins(dres, mres);
                    ratio_results[period][rkey] = rres;
                    const size_t nb = std::min({dres.size(), mres.size(), rres.size()});
                    for (size_t ib = 0; ib < nb; ++ib) {
                        ratio_csv << period << "," << rkey << "," << v.key << "," << ib << ","
                                  << std::setprecision(12)
                                  << rres[ib].x_low << "," << rres[ib].x_high << "," << rres[ib].x_center << ","
                                  << dres[ib].factor << "," << dres[ib].factor_err << ","
                                  << mres[ib].factor << "," << mres[ib].factor_err << ","
                                  << rres[ib].factor << "," << rres[ib].factor_err << ","
                                  << dres[ib].n_current_points << "," << mres[ib].n_current_points << "\n";
                    }
                }
            }
        }

        stored_data_results[v.key] = data_results;
        if (process_mc) {
            stored_mc_results[v.key] = mc_results;
            stored_ratio_results[v.key] = ratio_results;
        }

        draw_region_theta_canvas(v, "data", "Current-efficiency factor", data_results, 0.35, 1.45);
        if (process_mc) {
            draw_region_theta_canvas(v, "mc", "Current-efficiency factor", mc_results, 0.35, 1.45);
            draw_region_theta_canvas(v, "data_over_mc", "DATA/MC current-efficiency ratio",
                                     ratio_results, 0.50, 1.50);
        }

        // Model-selection/closure comparison on the relative current slopes.
        // M0 = one period-wide slope.
        // M1 = one independent baseline slope per photon region.
        // M2 = M1 plus one common linear dependence on this kinematic variable.
        for (const std::string& period : periods) {
            std::map<std::string, std::vector<DataAgg>> dmap;
            for (const PhotonRegionSpec& spec : PHOTON_REGION_ORDER)
                dmap[spec.key] = merged[period][spec.key][v.key];

            const auto dpts =
                relative_slope_points_from_data_theta_bins(dmap, v);
            const RelativeSlopeModelComparison dcmp =
                compare_relative_slope_models(dpts);
            if (dcmp.valid) {
                model_rows.push_back({"data", period, v.key, dcmp});
                model_csv << "data," << period << "," << v.key << ","
                          << dcmp.npoints << ","
                          << dcmp.chi2_m0 << "," << dcmp.ndf_m0 << ","
                          << dcmp.aic_m0 << "," << dcmp.bic_m0 << ","
                          << dcmp.chi2_m1 << "," << dcmp.ndf_m1 << ","
                          << dcmp.aic_m1 << "," << dcmp.bic_m1 << ","
                          << dcmp.chi2_m2 << "," << dcmp.ndf_m2 << ","
                          << dcmp.aic_m2 << "," << dcmp.bic_m2 << ","
                          << dcmp.gradient << "," << dcmp.gradient_err << ","
                          << dcmp.gradient_significance << ","
                          << (dcmp.chi2_m1 - dcmp.chi2_m2) << ","
                          << (dcmp.aic_m1 - dcmp.aic_m2) << ","
                          << (dcmp.bic_m1 - dcmp.bic_m2) << "\n";
            }

            if (process_mc) {
                std::map<std::string, std::vector<McKinematicBinAgg>> mmap;
                for (const PhotonRegionSpec& spec : PHOTON_REGION_ORDER)
                    mmap[spec.key] = merged_mc[period][spec.key][v.key];

                const auto mpts =
                    relative_slope_points_from_mc_theta_bins(mmap, v);
                const RelativeSlopeModelComparison mcmp =
                    compare_relative_slope_models(mpts);
                if (mcmp.valid) {
                    model_rows.push_back({"mc", period, v.key, mcmp});
                    model_csv << "mc," << period << "," << v.key << ","
                              << mcmp.npoints << ","
                              << mcmp.chi2_m0 << "," << mcmp.ndf_m0 << ","
                              << mcmp.aic_m0 << "," << mcmp.bic_m0 << ","
                              << mcmp.chi2_m1 << "," << mcmp.ndf_m1 << ","
                              << mcmp.aic_m1 << "," << mcmp.bic_m1 << ","
                              << mcmp.chi2_m2 << "," << mcmp.ndf_m2 << ","
                              << mcmp.aic_m2 << "," << mcmp.bic_m2 << ","
                              << mcmp.gradient << "," << mcmp.gradient_err << ","
                              << mcmp.gradient_significance << ","
                              << (mcmp.chi2_m1 - mcmp.chi2_m2) << ","
                              << (mcmp.aic_m1 - mcmp.aic_m2) << ","
                              << (mcmp.bic_m1 - mcmp.bic_m2) << "\n";
                }
            }
        }

        // Retain the dedicated e_theta pooled summary because it is especially
        // useful for the old Sp18-Out tracking-angle observation.
        if (v.key == "e_theta") {
            std::ofstream pooled_csv(odir + "/pooled_e_theta_relative_slope_model.csv");
            pooled_csv << "sample,period,theta_gradient_per_nA_per_deg,theta_gradient_stat,theta_gradient_percent_per_nA_per_deg,theta_gradient_percent_stat,chi2,ndf,npoints";
            for (const PhotonRegionSpec& spec : PHOTON_REGION_ORDER) {
                pooled_csv << "," << spec.key << "_theta_center_deg"
                           << "," << spec.key << "_relative_slope_per_nA"
                           << "," << spec.key << "_relative_slope_stat";
            }
            pooled_csv << "\n";

            std::vector<double> data_grad(periods.size(), std::numeric_limits<double>::quiet_NaN());
            std::vector<double> data_grad_err(periods.size(), std::numeric_limits<double>::quiet_NaN());
            std::vector<double> mc_grad(periods.size(), std::numeric_limits<double>::quiet_NaN());
            std::vector<double> mc_grad_err(periods.size(), std::numeric_limits<double>::quiet_NaN());

            for (size_t ip = 0; ip < periods.size(); ++ip) {
                const std::string& period = periods[ip];
                std::map<std::string, std::vector<DataAgg>> dmap;
                for (const PhotonRegionSpec& spec : PHOTON_REGION_ORDER)
                    dmap[spec.key] = merged[period][spec.key][v.key];

                const auto dpts = relative_slope_points_from_data_theta_bins(dmap, v);
                const auto dfit = fit_pooled_region_theta_relative_slopes(dpts);
                if (dfit.valid) {
                    production_e_theta_models[period] = dfit;
                    data_grad[ip] = 100.0 * dfit.theta_gradient;
                    data_grad_err[ip] = 100.0 * dfit.theta_gradient_err;
                }
                pooled_csv << "data," << period << ","
                           << dfit.theta_gradient << "," << dfit.theta_gradient_err << ","
                           << (100.0 * dfit.theta_gradient) << ","
                           << (100.0 * dfit.theta_gradient_err) << ","
                           << dfit.chi2 << "," << dfit.ndf << "," << dfit.npoints;
                for (int ir = 0; ir < 7; ++ir) {
                    pooled_csv << "," << dfit.region_centers[ir]
                               << "," << dfit.region_intercepts[ir]
                               << "," << dfit.region_intercept_err[ir];
                }
                pooled_csv << "\n";

                std::cout << "[current_dependence] pooled e_theta DATA relative-slope model "
                          << cfg.csv_channel << " " << period
                          << ": d(s_rel)/dtheta_e="
                          << 100.0 * dfit.theta_gradient << " +/- "
                          << 100.0 * dfit.theta_gradient_err
                          << " %/(nA deg), chi2/ndf="
                          << dfit.chi2 << "/" << dfit.ndf
                          << ", npoints=" << dfit.npoints << std::endl;

                if (process_mc) {
                    std::map<std::string, std::vector<McKinematicBinAgg>> mmap;
                    for (const PhotonRegionSpec& spec : PHOTON_REGION_ORDER)
                        mmap[spec.key] = merged_mc[period][spec.key][v.key];

                    const auto mpts = relative_slope_points_from_mc_theta_bins(mmap, v);
                    const auto mfit = fit_pooled_region_theta_relative_slopes(mpts);
                    if (mfit.valid) {
                        mc_grad[ip] = 100.0 * mfit.theta_gradient;
                        mc_grad_err[ip] = 100.0 * mfit.theta_gradient_err;
                    }
                    pooled_csv << "mc," << period << ","
                               << mfit.theta_gradient << "," << mfit.theta_gradient_err << ","
                               << (100.0 * mfit.theta_gradient) << ","
                               << (100.0 * mfit.theta_gradient_err) << ","
                               << mfit.chi2 << "," << mfit.ndf << "," << mfit.npoints;
                    for (int ir = 0; ir < 7; ++ir) {
                        pooled_csv << "," << mfit.region_centers[ir]
                                   << "," << mfit.region_intercepts[ir]
                                   << "," << mfit.region_intercept_err[ir];
                    }
                    pooled_csv << "\n";

                    const double diff = dfit.theta_gradient - mfit.theta_gradient;
                    const double diff_err = std::sqrt(
                        dfit.theta_gradient_err * dfit.theta_gradient_err +
                        mfit.theta_gradient_err * mfit.theta_gradient_err);
                    pooled_csv << "data_minus_mc," << period << ","
                               << diff << "," << diff_err << ","
                               << 100.0 * diff << "," << 100.0 * diff_err
                               << ",nan,0,0";
                    for (int ir = 0; ir < 7; ++ir) pooled_csv << ",nan,nan,nan";
                    pooled_csv << "\n";

                    std::cout << "[current_dependence] pooled e_theta DATA-MC gradient difference "
                              << cfg.csv_channel << " " << period
                              << ": " << 100.0 * diff << " +/- "
                              << 100.0 * diff_err << " %/(nA deg)" << std::endl;
                }
            }

            // Compact analysis-note summary of the pooled common theta_e
            // gradient. Units are percent change in the relative current slope
            // per nA per degree.
            TCanvas cg(("c_pooled_e_theta_gradient_" + cfg.output_token).c_str(), "", 1100, 750);
            cg.SetGridy();
            cg.SetLeftMargin(0.14);
            cg.SetBottomMargin(0.16);
            auto frame = std::make_unique<TH1D>(
                ("h_pooled_e_theta_gradient_" + cfg.output_token).c_str(),
                ";Run period;d s_{rel}/d#theta_{e}  (%/(nA deg))",
                (int)periods.size(), 0.5, periods.size() + 0.5);
            frame->SetDirectory(nullptr);
            for (size_t ip = 0; ip < periods.size(); ++ip)
                frame->GetXaxis()->SetBinLabel((int)ip + 1, periods[ip].c_str());

            double gymin = 0.0, gymax = 0.0;
            bool have_g = false;
            for (size_t ip = 0; ip < periods.size(); ++ip) {
                if (std::isfinite(data_grad[ip]) && std::isfinite(data_grad_err[ip])) {
                    gymin = have_g ? std::min(gymin, data_grad[ip] - data_grad_err[ip])
                                   : data_grad[ip] - data_grad_err[ip];
                    gymax = have_g ? std::max(gymax, data_grad[ip] + data_grad_err[ip])
                                   : data_grad[ip] + data_grad_err[ip];
                    have_g = true;
                }
                if (process_mc && std::isfinite(mc_grad[ip]) && std::isfinite(mc_grad_err[ip])) {
                    gymin = have_g ? std::min(gymin, mc_grad[ip] - mc_grad_err[ip])
                                   : mc_grad[ip] - mc_grad_err[ip];
                    gymax = have_g ? std::max(gymax, mc_grad[ip] + mc_grad_err[ip])
                                   : mc_grad[ip] + mc_grad_err[ip];
                    have_g = true;
                }
            }
            if (!have_g) { gymin = -0.05; gymax = 0.05; }
            const double gspan = std::max(0.02, gymax - gymin);
            frame->SetMinimum(gymin - 0.20 * gspan);
            frame->SetMaximum(gymax + 0.20 * gspan);
            frame->SetStats(0);
            frame->Draw();

            auto zero = std::make_unique<TLine>(0.5, 0.0, periods.size() + 0.5, 0.0);
            zero->SetLineStyle(2);
            zero->Draw("SAME");

            auto gd = std::make_unique<TGraphErrors>();
            auto gm = std::make_unique<TGraphErrors>();
            int nd = 0, nm = 0;
            for (size_t ip = 0; ip < periods.size(); ++ip) {
                if (std::isfinite(data_grad[ip]) && std::isfinite(data_grad_err[ip])) {
                    gd->SetPoint(nd, (double)ip + 1.0 - 0.08, data_grad[ip]);
                    gd->SetPointError(nd, 0.0, data_grad_err[ip]);
                    ++nd;
                }
                if (process_mc && std::isfinite(mc_grad[ip]) && std::isfinite(mc_grad_err[ip])) {
                    gm->SetPoint(nm, (double)ip + 1.0 + 0.08, mc_grad[ip]);
                    gm->SetPointError(nm, 0.0, mc_grad_err[ip]);
                    ++nm;
                }
            }
            gd->SetMarkerStyle(20);
            gd->SetMarkerSize(1.2);
            gd->SetLineWidth(2);
            gd->Draw("P SAME");
            if (process_mc) {
                gm->SetMarkerStyle(24);
                gm->SetMarkerSize(1.2);
                gm->SetLineWidth(2);
                gm->Draw("P SAME");
            }

            auto leg = std::make_unique<TLegend>(0.16, 0.78, 0.43, 0.90);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->AddEntry(gd.get(), "DATA", "pe");
            if (process_mc) leg->AddEntry(gm.get(), "MC", "pe");
            leg->Draw();

            cg.Modified();
            cg.Update();
            const std::string grad_path = odir + "/pooled_e_theta_gradient_data_mc.png";
            cg.SaveAs(grad_path.c_str());
            gSystem->CopyFile(
                grad_path.c_str(),
                (note_dir + "/" + cfg.output_token + "_pooled_e_theta_gradient_data_mc.png").c_str(),
                true);
        }
    }


    // ---------------------------------------------------------------------
    // FT-versus-FD relative current-response modeling in matched polar-angle
    // bins.  This is a current-efficiency diagnostic, not an absolute photon
    // reconstruction efficiency measurement.
    //
    // For each polar-angle bin:
    //   D = (DATA/MC)_FT / <(DATA/MC)_S1...S6>
    // where the FD denominator is an inverse-variance weighted mean of the six
    // sector ratios in the same angular bin.  Thus one polar-angle phase space
    // is held fixed at a time.  D=1 means FT and FD have the same relative
    // DATA/MC current response in that angular slice.
    // ---------------------------------------------------------------------
    std::ofstream ftfd_angle_csv(
        odir + "/ft_fd_data_over_mc_double_ratio_vs_polar_angles.csv");
    ftfd_angle_csv
        << "period,variable,bin,x_low,x_high,x_center,"
        << "ft_data_over_mc,ft_data_over_mc_stat,"
        << "fd_weighted_data_over_mc,fd_weighted_data_over_mc_stat,"
        << "ft_over_fd_double_ratio,ft_over_fd_double_ratio_stat,n_fd_sectors\n";

    for (const std::string& variable :
         {std::string("e_theta"), std::string("p_theta"), std::string("g_theta")}) {

        const KinematicVarConfig* vcfg = nullptr;
        for (const KinematicVarConfig& v : vars)
            if (v.key == variable) { vcfg = &v; break; }
        if (!vcfg) continue;

        auto ivar = stored_ratio_results.find(variable);
        if (ivar == stored_ratio_results.end()) continue;

        TCanvas cftfd(
            ("c_ft_fd_matched_" + cfg.output_token + "_" + variable).c_str(),
            "", 1800, 1200);
        cftfd.Divide(2, 2, 0.002, 0.002);
        std::vector<std::unique_ptr<TObject>> owned;

        for (size_t iperiod = 0; iperiod < periods.size(); ++iperiod) {
            const std::string& period = periods[iperiod];
            cftfd.cd((int)iperiod + 1);
            configure_production_pad((TPad*)gPad, 0.145, 0.035, 0.145, 0.075);

            auto frame = std::make_unique<TH1D>(
                ("h_ft_fd_matched_" + cfg.output_token + "_" + variable + "_" +
                 std::to_string(iperiod)).c_str(),
                (";" + vcfg->x_label +
                 ";(DATA/MC)_{FT} / (DATA/MC)_{FD}").c_str(),
                100, vcfg->edges.front(), vcfg->edges.back());
            frame->SetDirectory(nullptr);
            frame->SetMinimum(0.50);
            frame->SetMaximum(1.50);
            style_production_frame(frame.get(), 0.050, 0.050, 0.041);
            frame->Draw();
            owned.emplace_back(std::move(frame));

            auto unity = std::make_unique<TLine>(
                vcfg->edges.front(), 1.0, vcfg->edges.back(), 1.0);
            unity->SetLineStyle(2);
            unity->SetLineWidth(2);
            unity->SetLineColor(kGray + 1);
            unity->Draw("SAME");
            owned.emplace_back(std::move(unity));

            auto ip = ivar->second.find(period);
            if (ip != ivar->second.end()) {
                auto ift = ip->second.find("FT");
                if (ift != ip->second.end()) {
                    auto graph = std::make_unique<TGraphErrors>();
                    graph->SetMarkerStyle(20);
                    graph->SetMarkerSize(1.2);
                    graph->SetLineWidth(2);
                    int ng = 0;

                    const size_t nb = std::min(
                        ift->second.size(),
                        vcfg->edges.size() > 0 ? vcfg->edges.size() - 1 : size_t(0));

                    for (size_t ib = 0; ib < nb; ++ib) {
                        const KinematicBinResult& ft = ift->second[ib];
                        if (!(std::isfinite(ft.factor) && ft.factor > 0.0 &&
                              std::isfinite(ft.factor_err) && ft.factor_err >= 0.0)) {
                            continue;
                        }

                        double sw = 0.0;
                        double swy = 0.0;
                        int nfd = 0;
                        for (int is = 1; is <= 6; ++is) {
                            const std::string skey = "S" + std::to_string(is);
                            auto isec = ip->second.find(skey);
                            if (isec == ip->second.end() || ib >= isec->second.size()) continue;
                            const KinematicBinResult& sec = isec->second[ib];
                            if (!(std::isfinite(sec.factor) && sec.factor > 0.0 &&
                                  std::isfinite(sec.factor_err) && sec.factor_err > 0.0)) {
                                continue;
                            }
                            const double w = 1.0 / (sec.factor_err * sec.factor_err);
                            sw += w;
                            swy += w * sec.factor;
                            ++nfd;
                        }
                        if (!(sw > 0.0) || nfd == 0) continue;

                        const double fd = swy / sw;
                        const double sfd = std::sqrt(1.0 / sw);
                        double sdouble = std::numeric_limits<double>::quiet_NaN();
                        const double d = ratio_with_error(
                            ft.factor, ft.factor_err, fd, sfd, sdouble);
                        if (!std::isfinite(d) || !std::isfinite(sdouble)) continue;

                        graph->SetPoint(ng, ft.x_center, d);
                        graph->SetPointError(ng, 0.0, sdouble);
                        ++ng;

                        ftfd_angle_csv
                            << period << "," << variable << "," << ib << ","
                            << ft.x_low << "," << ft.x_high << "," << ft.x_center << ","
                            << ft.factor << "," << ft.factor_err << ","
                            << fd << "," << sfd << ","
                            << d << "," << sdouble << "," << nfd << "\n";
                    }

                    graph->Draw("PE SAME");
                    owned.emplace_back(std::move(graph));
                }
            }

            owned.emplace_back(
                make_period_label(period, 0.94, 0.86, 33, 0.052));
        }

        cftfd.cd(0);
        auto title = std::make_unique<TLatex>();
        title->SetNDC();
        title->SetTextAlign(22);
        title->SetTextFont(62);
        title->SetTextSize(0.029);
        title->DrawLatex(
            0.5, 0.975,
            (cfg.title + ": FT/FD relative DATA/MC current response at matched " +
             vcfg->title).c_str());
        owned.emplace_back(std::move(title));

        cftfd.Modified();
        cftfd.Update();
        const std::string path =
            odir + "/ft_fd_data_over_mc_double_ratio_vs_" + variable + ".png";
        cftfd.SaveAs(path.c_str());
        gSystem->CopyFile(
            path.c_str(),
            (note_dir + "/" + cfg.output_token +
             "_ft_fd_data_over_mc_double_ratio_vs_" + variable + ".png").c_str(),
            true);
    }

    // ---------------------------------------------------------------------
    // Compact model decision table.
    //
    // Keep the model space deliberately small:
    //   M1               : regional-only response
    //   theta_e linear   : M1 + one common linear theta_e gradient
    //   theta_e quadratic: M1 + common centered linear and quadratic theta_e terms
    //   best polar angle : best one-variable linear extension among
    //                      theta_e, theta_p, and theta_gamma only
    //
    // This is intended to answer whether modest added complexity actually buys
    // useful closure, not to turn the current-response fit into a large
    // multidimensional regression.
    // ---------------------------------------------------------------------
    std::ofstream compact_csv(odir + "/compact_current_response_model_comparison.csv");
    compact_csv
        << "sample,period,regional_chi2,regional_ndf,regional_bic,"
        << "theta_e_linear_chi2,theta_e_linear_ndf,theta_e_linear_bic,theta_e_linear_delta_bic,"
        << "theta_e_quadratic_chi2,theta_e_quadratic_ndf,theta_e_quadratic_bic,theta_e_quadratic_delta_bic,"
        << "theta_e_quadratic_linear,theta_e_quadratic_linear_stat,"
        << "theta_e_quadratic_coefficient,theta_e_quadratic_coefficient_stat,"
        << "best_polar_angle_variable,best_polar_angle_chi2,best_polar_angle_ndf,best_polar_angle_bic,best_polar_angle_delta_bic\n";

    auto find_var_config = [&](const std::string& key) -> const KinematicVarConfig* {
        for (const KinematicVarConfig& v : vars)
            if (v.key == key) return &v;
        return nullptr;
    };

    auto find_model_row =
        [&](const std::string& sample,
            const std::string& period,
            const std::string& variable) -> const CandidateModelRow* {
        for (const CandidateModelRow& row : model_rows) {
            if (row.sample == sample && row.period == period &&
                row.variable == variable && row.cmp.valid) return &row;
        }
        return nullptr;
    };

    auto best_polar_angle_row =
        [&](const std::string& sample,
            const std::string& period) -> const CandidateModelRow* {
        const CandidateModelRow* best = nullptr;
        double best_delta_bic = -std::numeric_limits<double>::infinity();
        for (const CandidateModelRow& row : model_rows) {
            if (row.sample != sample || row.period != period || !row.cmp.valid) continue;
            if (!is_polar_angle_model_variable(row.variable)) continue;
            const double dbic = row.cmp.bic_m1 - row.cmp.bic_m2;
            if (!best || dbic > best_delta_bic) {
                best = &row;
                best_delta_bic = dbic;
            }
        }
        return best;
    };

    const KinematicVarConfig* theta_e_cfg = find_var_config("e_theta");

    struct CompactDecisionRow {
        std::string sample;
        std::string period;
        double regional_bic = std::numeric_limits<double>::quiet_NaN();
        double theta_linear_bic = std::numeric_limits<double>::quiet_NaN();
        double theta_quadratic_bic = std::numeric_limits<double>::quiet_NaN();
        std::string best_polar_angle_variable;
        double best_polar_angle_delta_bic = std::numeric_limits<double>::quiet_NaN();
    };
    std::vector<CompactDecisionRow> compact_rows;

    auto write_compact_row =
        [&](const std::string& sample,
            const std::string& period) {

        const CandidateModelRow* theta_linear =
            find_model_row(sample, period, "e_theta");
        const CandidateModelRow* best =
            best_polar_angle_row(sample, period);

        PooledRegionQuadraticFit theta_quad;
        if (theta_e_cfg) {
            if (sample == "data") {
                std::map<std::string, std::vector<DataAgg>> dmap;
                for (const PhotonRegionSpec& spec : PHOTON_REGION_ORDER)
                    dmap[spec.key] = merged[period][spec.key][theta_e_cfg->key];
                theta_quad = fit_pooled_region_quadratic_relative_slopes(
                    relative_slope_points_from_data_theta_bins(dmap, *theta_e_cfg));
            } else if (sample == "mc" && process_mc) {
                std::map<std::string, std::vector<McKinematicBinAgg>> mmap;
                for (const PhotonRegionSpec& spec : PHOTON_REGION_ORDER)
                    mmap[spec.key] = merged_mc[period][spec.key][theta_e_cfg->key];
                theta_quad = fit_pooled_region_quadratic_relative_slopes(
                    relative_slope_points_from_mc_theta_bins(mmap, *theta_e_cfg));
            }
        }

        const double m1_chi2 = theta_linear ? theta_linear->cmp.chi2_m1
                                            : std::numeric_limits<double>::quiet_NaN();
        const int m1_ndf = theta_linear ? theta_linear->cmp.ndf_m1 : 0;
        const double m1_bic = theta_linear ? theta_linear->cmp.bic_m1
                                           : std::numeric_limits<double>::quiet_NaN();

        CompactDecisionRow compact_row;
        compact_row.sample = sample;
        compact_row.period = period;
        compact_row.regional_bic = m1_bic;
        compact_row.theta_linear_bic =
            theta_linear ? theta_linear->cmp.bic_m2
                         : std::numeric_limits<double>::quiet_NaN();
        compact_row.theta_quadratic_bic =
            theta_quad.valid ? theta_quad.bic
                             : std::numeric_limits<double>::quiet_NaN();
        compact_row.best_polar_angle_variable = best ? best->variable : "";
        compact_row.best_polar_angle_delta_bic =
            best ? (best->cmp.bic_m1 - best->cmp.bic_m2)
                 : std::numeric_limits<double>::quiet_NaN();
        compact_rows.push_back(compact_row);

        if (sample == "data") {
            std::cout << "[current_dependence] compact model comparison "
                      << cfg.csv_channel << " " << period
                      << ": regional chi2/ndf=" << m1_chi2 << "/" << m1_ndf
                      << ", theta_e linear DeltaBIC="
                      << (theta_linear
                              ? theta_linear->cmp.bic_m1 - theta_linear->cmp.bic_m2
                              : std::numeric_limits<double>::quiet_NaN())
                      << ", theta_e quadratic DeltaBIC="
                      << ((theta_quad.valid && std::isfinite(m1_bic))
                              ? m1_bic - theta_quad.bic
                              : std::numeric_limits<double>::quiet_NaN())
                      << ", best polar angle="
                      << (best ? best->variable : std::string("n/a"))
                      << " (DeltaBIC="
                      << (best ? best->cmp.bic_m1 - best->cmp.bic_m2
                               : std::numeric_limits<double>::quiet_NaN())
                      << ")"
                      << std::endl;
        }

        compact_csv << sample << "," << period << ","
                    << m1_chi2 << "," << m1_ndf << "," << m1_bic << ","
                    << (theta_linear ? theta_linear->cmp.chi2_m2 : std::numeric_limits<double>::quiet_NaN()) << ","
                    << (theta_linear ? theta_linear->cmp.ndf_m2 : 0) << ","
                    << (theta_linear ? theta_linear->cmp.bic_m2 : std::numeric_limits<double>::quiet_NaN()) << ","
                    << (theta_linear ? theta_linear->cmp.bic_m1 - theta_linear->cmp.bic_m2
                                     : std::numeric_limits<double>::quiet_NaN()) << ","
                    << theta_quad.chi2 << "," << theta_quad.ndf << "," << theta_quad.bic << ","
                    << ((theta_quad.valid && std::isfinite(m1_bic)) ? m1_bic - theta_quad.bic
                                                                  : std::numeric_limits<double>::quiet_NaN()) << ","
                    << theta_quad.linear << "," << theta_quad.linear_err << ","
                    << theta_quad.quadratic << "," << theta_quad.quadratic_err << ","
                    << (best ? best->variable : "") << ","
                    << (best ? best->cmp.chi2_m2 : std::numeric_limits<double>::quiet_NaN()) << ","
                    << (best ? best->cmp.ndf_m2 : 0) << ","
                    << (best ? best->cmp.bic_m2 : std::numeric_limits<double>::quiet_NaN()) << ","
                    << (best ? best->cmp.bic_m1 - best->cmp.bic_m2
                             : std::numeric_limits<double>::quiet_NaN())
                    << "\n";
    };

    for (const std::string& period : periods) {
        write_compact_row("data", period);
        if (process_mc) write_compact_row("mc", period);
    }

    // Compact production-quality model decision plot for DATA.  Values are
    // DeltaBIC relative to the regional-only model; larger positive values
    // indicate that the modest extension is preferred.
    TCanvas ccompact(("c_compact_model_choice_" + cfg.output_token).c_str(),
                     "", 1800, 1200);
    ccompact.Divide(2, 2, 0.002, 0.002);
    std::vector<std::unique_ptr<TObject>> compact_owned;

    for (size_t ip = 0; ip < periods.size(); ++ip) {
        ccompact.cd((int)ip + 1);
        configure_production_pad((TPad*)gPad, 0.145, 0.035, 0.17, 0.075);

        const CompactDecisionRow* row = nullptr;
        for (const CompactDecisionRow& r : compact_rows) {
            if (r.sample == "data" && r.period == periods[ip]) {
                row = &r;
                break;
            }
        }

        std::vector<std::string> labels = {
            "#theta_{e} linear",
            "#theta_{e} quadratic",
            row && !row->best_polar_angle_variable.empty()
                ? ("best polar angle: " + pretty_kinematic_key(row->best_polar_angle_variable))
                : "best polar angle"
        };
        std::vector<double> values(3, std::numeric_limits<double>::quiet_NaN());
        if (row && std::isfinite(row->regional_bic)) {
            if (std::isfinite(row->theta_linear_bic))
                values[0] = row->regional_bic - row->theta_linear_bic;
            if (std::isfinite(row->theta_quadratic_bic))
                values[1] = row->regional_bic - row->theta_quadratic_bic;
            if (std::isfinite(row->best_polar_angle_delta_bic))
                values[2] = row->best_polar_angle_delta_bic;
        }

        double ymin = -5.0;
        double ymax = 5.0;
        for (double v : values) {
            if (!std::isfinite(v)) continue;
            ymin = std::min(ymin, v);
            ymax = std::max(ymax, v);
        }
        const double span = std::max(5.0, ymax - ymin);

        auto frame = std::make_unique<TH1D>(
            ("h_compact_model_choice_" + cfg.output_token + "_" + std::to_string(ip)).c_str(),
            ";Candidate model;#DeltaBIC relative to regional-only",
            3, 0.5, 3.5);
        frame->SetDirectory(nullptr);
        for (int ib = 0; ib < 3; ++ib)
            frame->GetXaxis()->SetBinLabel(ib + 1, labels[ib].c_str());
        frame->SetMinimum(std::min(-5.0, ymin - 0.12 * span));
        frame->SetMaximum(std::max(5.0, ymax + 0.15 * span));
        style_production_frame(frame.get(), 0.050, 0.050, 0.039);
        frame->GetXaxis()->SetLabelSize(0.035);
        frame->Draw();
        compact_owned.emplace_back(std::move(frame));

        auto zero = std::make_unique<TLine>(0.5, 0.0, 3.5, 0.0);
        zero->SetLineStyle(2);
        zero->SetLineWidth(2);
        zero->SetLineColor(kGray + 1);
        zero->Draw("SAME");
        compact_owned.emplace_back(std::move(zero));

        auto graph = std::make_unique<TGraph>();
        int ng = 0;
        for (int iv = 0; iv < 3; ++iv) {
            if (std::isfinite(values[iv]))
                graph->SetPoint(ng++, iv + 1.0, values[iv]);
        }
        graph->SetMarkerStyle(20);
        graph->SetMarkerSize(1.4);
        graph->SetLineWidth(2);
        graph->Draw("P SAME");
        compact_owned.emplace_back(std::move(graph));

        compact_owned.emplace_back(
            make_period_label(periods[ip], 0.94, 0.86, 33, 0.052));
    }

    ccompact.cd(0);
    auto compact_title = std::make_unique<TLatex>();
    compact_title->SetNDC();
    compact_title->SetTextAlign(22);
    compact_title->SetTextFont(62);
    compact_title->SetTextSize(0.030);
    compact_title->DrawLatex(
        0.5, 0.972,
        (cfg.title + ": compact polar-angle current-response model comparison").c_str());
    compact_owned.emplace_back(std::move(compact_title));

    ccompact.Modified();
    ccompact.Update();
    const std::string compact_plot =
        odir + "/compact_current_response_model_comparison.png";
    ccompact.SaveAs(compact_plot.c_str());
    gSystem->CopyFile(
        compact_plot.c_str(),
        (note_dir + "/" + cfg.output_token +
         "_compact_current_response_model_comparison.png").c_str(),
        true);

    // ---------------------------------------------------------------------
    // Targeted production-quality evidence canvases.
    // These put DATA, reconstructed MC, and DATA/MC next to one another for the
    // two effects that presently matter most:
    //   Sp18 Out : theta_e
    //   Fa18 Out : theta_p
    // ---------------------------------------------------------------------
    auto draw_targeted_evidence =
        [&](const std::string& period,
            const std::string& variable,
            const std::string& output_token,
            const std::string& headline) {

        auto idv = stored_data_results.find(variable);
        auto imv = stored_mc_results.find(variable);
        auto irv = stored_ratio_results.find(variable);
        const KinematicVarConfig* vcfg = find_var_config(variable);
        if (!vcfg || idv == stored_data_results.end() ||
            imv == stored_mc_results.end() || irv == stored_ratio_results.end()) return;

        auto idp = idv->second.find(period);
        auto imp = imv->second.find(period);
        auto irp = irv->second.find(period);
        if (idp == idv->second.end() || imp == imv->second.end() || irp == irv->second.end()) return;

        struct PanelSpec {
            std::string label;
            const std::map<std::string, std::vector<KinematicBinResult>>* results = nullptr;
            bool unity = false;
        };
        std::vector<PanelSpec> panels = {
            {"DATA", &idp->second, false},
            {"MC", &imp->second, false},
            {"DATA/MC", &irp->second, true}
        };

        TCanvas c(("c_targeted_" + cfg.output_token + "_" + output_token).c_str(),
                  "", 2100, 720);
        c.Divide(3, 1, 0.002, 0.002);
        std::vector<std::unique_ptr<TObject>> owned;

        for (int ipanel = 0; ipanel < 3; ++ipanel) {
            c.cd(ipanel + 1);
            configure_production_pad((TPad*)gPad, 0.15, 0.035, 0.15, 0.085);

            const auto yrange = production_efficiency_plot_range(
                panels[ipanel].unity ? "data_over_mc" : "factor");
            const double lo = yrange.first;
            const double hi = yrange.second;

            auto frame = std::make_unique<TH1D>(
                ("h_targeted_" + cfg.output_token + "_" + output_token + "_" +
                 std::to_string(ipanel)).c_str(),
                (";" + vcfg->x_label + ";" +
                 (panels[ipanel].unity ? "DATA/MC current-efficiency ratio"
                                      : "Current-efficiency factor")).c_str(),
                100, vcfg->edges.front(), vcfg->edges.back());
            frame->SetDirectory(nullptr);
            frame->SetMinimum(lo);
            frame->SetMaximum(hi);
            style_production_frame(frame.get(), 0.050, 0.050, 0.041);
            frame->Draw();
            owned.emplace_back(std::move(frame));

            if (panels[ipanel].unity) {
                auto line = std::make_unique<TLine>(
                    vcfg->edges.front(), 1.0, vcfg->edges.back(), 1.0);
                line->SetLineStyle(2);
                line->SetLineWidth(2);
                line->SetLineColor(kGray + 1);
                line->Draw("SAME");
                owned.emplace_back(std::move(line));
            }

            std::vector<TGraphErrors*> legend_graphs;
            for (int ir = 0; ir < 7; ++ir) {
                auto it = panels[ipanel].results->find(PHOTON_REGION_ORDER[ir].key);
                if (it == panels[ipanel].results->end()) continue;
                auto graph = std::unique_ptr<TGraphErrors>(
                    make_region_theta_factor_graph(it->second, photon_region_color(ir)));
                graph->SetMarkerStyle(photon_region_marker(ir));
                graph->SetMarkerSize(1.15);
                graph->SetLineWidth(2);
                graph->Draw("P SAME");
                if (ipanel == 0) legend_graphs.push_back(graph.get());
                owned.emplace_back(std::move(graph));
            }

            auto panel_label = std::make_unique<TLatex>();
            panel_label->SetNDC();
            panel_label->SetTextFont(62);
            panel_label->SetTextSize(0.052);
            panel_label->DrawLatex(0.17, 0.84, panels[ipanel].label.c_str());
            owned.emplace_back(std::move(panel_label));

            if (ipanel == 0) {
                auto leg = std::make_unique<TLegend>(0.17, 0.71, 0.80, 0.87);
                leg->SetNColumns(4);
                leg->SetBorderSize(0);
                leg->SetFillStyle(0);
                leg->SetTextFont(42);
                leg->SetTextSize(0.033);
                for (size_t ir = 0; ir < legend_graphs.size() && ir < PHOTON_REGION_ORDER.size(); ++ir)
                    leg->AddEntry(legend_graphs[ir], PHOTON_REGION_ORDER[ir].label.c_str(), "pe");
                leg->Draw();
                owned.emplace_back(std::move(leg));
            }
        }

        c.cd(0);
        auto title = std::make_unique<TLatex>();
        title->SetNDC();
        title->SetTextAlign(22);
        title->SetTextFont(62);
        title->SetTextSize(0.033);
        title->DrawLatex(0.5, 0.988, (period + ": " + headline).c_str());
        owned.emplace_back(std::move(title));

        c.Modified();
        c.Update();
        const std::string path = odir + "/" + output_token + ".png";
        c.SaveAs(path.c_str());
        gSystem->CopyFile(
            path.c_str(),
            (note_dir + "/" + cfg.output_token + "_" + output_token + ".png").c_str(),
            true);
    };

    if (process_mc && cfg.channel == Channel::DVCS) {
        draw_targeted_evidence(
            "Sp18 Out", "e_theta",
            "sp18_out_theta_e_current_response",
            "residual current response versus electron polar angle");
        draw_targeted_evidence(
            "Fa18 Out", "p_theta",
            "fa18_out_theta_p_current_response",
            "residual current response versus proton polar angle");
    }

    // ---------------------------------------------------------------------
    // Closure pulls for the two leading outbending effects.
    // Left: regional-only model (M1).
    // Right: regional + one common linear dependence (M2).
    // This makes it immediately visible whether the extra term removes a
    // coherent residual trend rather than merely winning an information
    // criterion numerically.
    // ---------------------------------------------------------------------
    std::ofstream closure_csv(odir + "/targeted_model_closure_pulls.csv");
    closure_csv
        << "period,variable,region,x,relative_slope,relative_slope_stat,"
        << "regional_prediction,linear_prediction,regional_pull,linear_pull\n";

    auto draw_model_closure =
        [&](const std::string& period,
            const std::string& variable,
            const std::string& output_token) {

        const KinematicVarConfig* vcfg = find_var_config(variable);
        if (!vcfg) return;

        std::map<std::string, std::vector<DataAgg>> dmap;
        for (const PhotonRegionSpec& spec : PHOTON_REGION_ORDER)
            dmap[spec.key] = merged[period][spec.key][variable];

        const std::vector<RelativeSlopePoint> pts =
            relative_slope_points_from_data_theta_bins(dmap, *vcfg);
        const PooledRegionThetaSlopeFit linear =
            fit_pooled_region_theta_relative_slopes(pts);
        if (!linear.valid) return;

        std::array<double, 7> sw{};
        std::array<double, 7> swy{};
        for (const RelativeSlopePoint& pt : pts) {
            if (pt.region_index < 0 || pt.region_index >= 7 ||
                !std::isfinite(pt.slope) || !std::isfinite(pt.slope_err) ||
                !(pt.slope_err > 0.0)) continue;
            const double w = 1.0 / (pt.slope_err * pt.slope_err);
            sw[pt.region_index] += w;
            swy[pt.region_index] += w * pt.slope;
        }
        std::array<double, 7> regional_mean{};
        for (int ir = 0; ir < 7; ++ir)
            regional_mean[ir] = sw[ir] > 0.0 ? swy[ir] / sw[ir] : 0.0;

        TCanvas c(("c_closure_" + cfg.output_token + "_" + output_token).c_str(),
                  "", 1500, 720);
        c.Divide(2, 1, 0.002, 0.002);
        std::vector<std::unique_ptr<TObject>> owned;

        for (const RelativeSlopePoint& pt : pts) {
            if (pt.region_index < 0 || pt.region_index >= 7 ||
                !std::isfinite(pt.slope_err) || !(pt.slope_err > 0.0)) continue;
            const double p1 = (pt.slope - regional_mean[pt.region_index]) / pt.slope_err;
            const double pred2 =
                linear.region_intercepts[pt.region_index] +
                linear.theta_gradient * (pt.x - linear.region_centers[pt.region_index]);
            const double p2 = (pt.slope - pred2) / pt.slope_err;
            closure_csv << period << "," << variable << ","
                        << PHOTON_REGION_ORDER[pt.region_index].key << ","
                        << pt.x << "," << pt.slope << "," << pt.slope_err << ","
                        << regional_mean[pt.region_index] << "," << pred2 << ","
                        << p1 << "," << p2 << "\n";
        }
        const double pull_lim = production_pull_limit();

        for (int imod = 0; imod < 2; ++imod) {
            c.cd(imod + 1);
            configure_production_pad((TPad*)gPad, 0.14, 0.035, 0.15, 0.085);

            auto frame = std::make_unique<TH1D>(
                ("h_closure_" + cfg.output_token + "_" + output_token + "_" +
                 std::to_string(imod)).c_str(),
                (";" + vcfg->x_label + ";Residual pull").c_str(),
                100, vcfg->edges.front(), vcfg->edges.back());
            frame->SetDirectory(nullptr);
            frame->SetMinimum(-pull_lim);
            frame->SetMaximum(pull_lim);
            style_production_frame(frame.get(), 0.050, 0.050, 0.041);
            frame->Draw();
            owned.emplace_back(std::move(frame));

            auto zero = std::make_unique<TLine>(
                vcfg->edges.front(), 0.0, vcfg->edges.back(), 0.0);
            zero->SetLineStyle(2);
            zero->SetLineWidth(2);
            zero->SetLineColor(kGray + 1);
            zero->Draw("SAME");
            owned.emplace_back(std::move(zero));

            std::array<std::unique_ptr<TGraph>, 7> graphs;
            for (int ir = 0; ir < 7; ++ir) {
                graphs[ir] = std::make_unique<TGraph>();
                graphs[ir]->SetMarkerStyle(photon_region_marker(ir));
                graphs[ir]->SetMarkerColor(photon_region_color(ir));
                graphs[ir]->SetLineColor(photon_region_color(ir));
                graphs[ir]->SetMarkerSize(1.15);
            }

            std::array<int, 7> ngraph{};
            for (const RelativeSlopePoint& pt : pts) {
                if (pt.region_index < 0 || pt.region_index >= 7 ||
                    !std::isfinite(pt.slope_err) || !(pt.slope_err > 0.0)) continue;
                double pull = 0.0;
                if (imod == 0) {
                    pull = (pt.slope - regional_mean[pt.region_index]) / pt.slope_err;
                } else {
                    const double pred =
                        linear.region_intercepts[pt.region_index] +
                        linear.theta_gradient *
                        (pt.x - linear.region_centers[pt.region_index]);
                    pull = (pt.slope - pred) / pt.slope_err;
                }
                graphs[pt.region_index]->SetPoint(
                    ngraph[pt.region_index]++, pt.x, pull);
            }

            std::vector<TGraph*> legend_graphs;
            for (int ir = 0; ir < 7; ++ir) {
                graphs[ir]->Draw("P SAME");
                if (imod == 0) legend_graphs.push_back(graphs[ir].get());
                owned.emplace_back(std::move(graphs[ir]));
            }

            auto model_label = std::make_unique<TLatex>();
            model_label->SetNDC();
            model_label->SetTextFont(62);
            model_label->SetTextSize(0.050);
            model_label->DrawLatex(
                0.17, 0.84,
                imod == 0 ? "Regional only (M1)" : "Regional + linear term (M2)");
            owned.emplace_back(std::move(model_label));

            if (imod == 0) {
                auto leg = std::make_unique<TLegend>(0.17, 0.73, 0.80, 0.87);
                leg->SetNColumns(4);
                leg->SetBorderSize(0);
                leg->SetFillStyle(0);
                leg->SetTextFont(42);
                leg->SetTextSize(0.032);
                for (size_t ir = 0; ir < legend_graphs.size(); ++ir)
                    leg->AddEntry(legend_graphs[ir], PHOTON_REGION_ORDER[ir].label.c_str(), "p");
                leg->Draw();
                owned.emplace_back(std::move(leg));
            }
        }

        c.cd(0);
        auto title = std::make_unique<TLatex>();
        title->SetNDC();
        title->SetTextAlign(22);
        title->SetTextFont(62);
        title->SetTextSize(0.033);
        title->DrawLatex(
            0.5, 0.972,
            (period + ": current-response closure versus " +
             pretty_kinematic_key(variable)).c_str());
        owned.emplace_back(std::move(title));

        c.Modified();
        c.Update();
        const std::string path = odir + "/" + output_token + ".png";
        c.SaveAs(path.c_str());
        gSystem->CopyFile(
            path.c_str(),
            (note_dir + "/" + cfg.output_token + "_" + output_token + ".png").c_str(),
            true);
    };

    if (cfg.channel == Channel::DVCS) {
        draw_model_closure(
            "Sp18 Out", "e_theta",
            "sp18_out_theta_e_model_closure");
        draw_model_closure(
            "Fa18 Out", "p_theta",
            "fa18_out_theta_p_model_closure");
    }

    // ---------------------------------------------------------------------
    // Explicit lowest-bin robustness check for the apparent Fa18 Out theta_p
    // preference.  This is diagnostic only: no low-statistics point is silently
    // removed from the ordinary plots or production response.  We simply repeat
    // the same nested M1/M2 comparison after omitting the configured first
    // theta_p interval (10--18 deg) from every photon region.
    // ---------------------------------------------------------------------
    if (cfg.channel == Channel::DVCS) {
        const KinematicVarConfig* pcfg = find_var_config("p_theta");
        if (pcfg && pcfg->edges.size() >= 3) {
            std::map<std::string, std::vector<DataAgg>> dmap;
            for (const PhotonRegionSpec& spec : PHOTON_REGION_ORDER)
                dmap[spec.key] = merged["Fa18 Out"][spec.key]["p_theta"];

            const std::vector<RelativeSlopePoint> full_pts =
                relative_slope_points_from_data_theta_bins(dmap, *pcfg);
            std::vector<RelativeSlopePoint> trimmed_pts;
            for (const RelativeSlopePoint& pt : full_pts) {
                if (pt.x >= pcfg->edges[1]) trimmed_pts.push_back(pt);
            }

            const RelativeSlopeModelComparison full_cmp =
                compare_relative_slope_models(full_pts);
            const RelativeSlopeModelComparison trim_cmp =
                compare_relative_slope_models(trimmed_pts);

            std::ofstream rob_csv(
                odir + "/fa18_out_theta_p_lowest_bin_robustness.csv");
            rob_csv
                << "selection,npoints,lowest_theta_deg,regional_chi2,regional_ndf,"
                << "linear_chi2,linear_ndf,delta_chi2,delta_bic,gradient,"
                << "gradient_stat,gradient_significance\n";
            auto write_rob = [&](const std::string& sel,
                                 const RelativeSlopeModelComparison& cmp,
                                 int np,
                                 double low) {
                rob_csv << sel << "," << np << "," << low << ","
                        << cmp.chi2_m1 << "," << cmp.ndf_m1 << ","
                        << cmp.chi2_m2 << "," << cmp.ndf_m2 << ","
                        << cmp.chi2_m1 - cmp.chi2_m2 << ","
                        << cmp.bic_m1 - cmp.bic_m2 << ","
                        << cmp.gradient << "," << cmp.gradient_err << ","
                        << cmp.gradient_significance << "\n";
            };
            write_rob("full", full_cmp, (int)full_pts.size(), pcfg->edges.front());
            write_rob("omit_10_18_deg", trim_cmp, (int)trimmed_pts.size(), pcfg->edges[1]);

            // Four-panel figure: the actual factors, the event statistics that
            // expose the pathological edge bin, the regional-only pulls, and
            // the full-vs-omit-first-bin BIC comparison.
            auto iv = stored_data_results.find("p_theta");
            if (iv != stored_data_results.end()) {
                auto iper = iv->second.find("Fa18 Out");
                if (iper != iv->second.end()) {
                    TCanvas crob(
                        ("c_fa18_out_theta_p_robustness_" + cfg.output_token).c_str(),
                        "", 1800, 1200);
                    crob.Divide(2, 2, 0.002, 0.002);
                    std::vector<std::unique_ptr<TObject>> rob_owned;

                    // Panel 1: current-efficiency factors, fixed production scale.
                    crob.cd(1);
                    configure_production_pad((TPad*)gPad, 0.145, 0.035, 0.145, 0.085);
                    auto f1 = std::make_unique<TH1D>(
                        ("h_rob_factor_" + cfg.output_token).c_str(),
                        ";#theta_{p} (deg);DATA current-efficiency factor",
                        100, pcfg->edges.front(), pcfg->edges.back());
                    f1->SetDirectory(nullptr);
                    f1->SetMinimum(0.20);
                    f1->SetMaximum(1.60);
                    style_production_frame(f1.get(), 0.050, 0.050, 0.041);
                    f1->Draw();
                    rob_owned.emplace_back(std::move(f1));

                    std::vector<TGraphErrors*> leggraphs;
                    for (int ir = 0; ir < 7; ++ir) {
                        auto ireg = iper->second.find(PHOTON_REGION_ORDER[ir].key);
                        if (ireg == iper->second.end()) continue;
                        auto g = std::unique_ptr<TGraphErrors>(
                            make_region_theta_factor_graph(
                                ireg->second, photon_region_color(ir)));
                        g->SetMarkerStyle(photon_region_marker(ir));
                        g->SetMarkerSize(1.15);
                        g->SetLineWidth(2);
                        g->Draw("P SAME");
                        leggraphs.push_back(g.get());
                        rob_owned.emplace_back(std::move(g));
                    }
                    auto l1 = std::make_unique<TLegend>(0.17, 0.70, 0.80, 0.84);
                    l1->SetNColumns(4);
                    l1->SetBorderSize(0);
                    l1->SetFillStyle(0);
                    l1->SetTextFont(42);
                    l1->SetTextSize(0.032);
                    for (size_t ir = 0; ir < leggraphs.size(); ++ir)
                        l1->AddEntry(
                            leggraphs[ir], PHOTON_REGION_ORDER[ir].label.c_str(), "pe");
                    l1->Draw();
                    rob_owned.emplace_back(std::move(l1));
                    auto p1lab = std::make_unique<TLatex>();
                    p1lab->SetNDC(); p1lab->SetTextFont(62); p1lab->SetTextSize(0.047);
                    p1lab->DrawLatex(0.17, 0.88, "Current-response factors");
                    rob_owned.emplace_back(std::move(p1lab));

                    // Panel 2: event statistics by angular bin.
                    crob.cd(2);
                    configure_production_pad((TPad*)gPad, 0.145, 0.035, 0.145, 0.085);
                    gPad->SetLogy();
                    auto f2 = std::make_unique<TH1D>(
                        ("h_rob_counts_" + cfg.output_token).c_str(),
                        ";#theta_{p} (deg);Selected DATA events",
                        100, pcfg->edges.front(), pcfg->edges.back());
                    f2->SetDirectory(nullptr);
                    f2->SetMinimum(0.7);
                    f2->SetMaximum(1.0e6);
                    style_production_frame(f2.get(), 0.050, 0.050, 0.041);
                    f2->Draw();
                    rob_owned.emplace_back(std::move(f2));
                    for (int ir = 0; ir < 7; ++ir) {
                        auto ireg = iper->second.find(PHOTON_REGION_ORDER[ir].key);
                        if (ireg == iper->second.end()) continue;
                        auto g = std::make_unique<TGraph>();
                        int n = 0;
                        for (const KinematicBinResult& b : ireg->second) {
                            if (!(b.total_counts > 0.0)) continue;
                            g->SetPoint(n++, b.x_center, b.total_counts);
                        }
                        g->SetMarkerStyle(photon_region_marker(ir));
                        g->SetMarkerColor(photon_region_color(ir));
                        g->SetLineColor(photon_region_color(ir));
                        g->SetMarkerSize(1.15);
                        g->Draw("P SAME");
                        rob_owned.emplace_back(std::move(g));
                    }
                    auto edge = std::make_unique<TLine>(
                        pcfg->edges[1], 0.7, pcfg->edges[1], 1.0e6);
                    edge->SetLineStyle(2); edge->SetLineWidth(2);
                    edge->SetLineColor(kGray + 1); edge->Draw("SAME");
                    rob_owned.emplace_back(std::move(edge));
                    auto p2lab = std::make_unique<TLatex>();
                    p2lab->SetNDC(); p2lab->SetTextFont(62); p2lab->SetTextSize(0.047);
                    p2lab->DrawLatex(0.17, 0.88, "Statistics by #theta_{p} bin");
                    rob_owned.emplace_back(std::move(p2lab));

                    // Panel 3: regional-only pulls, deliberately clipped at +/-6.
                    crob.cd(3);
                    configure_production_pad((TPad*)gPad, 0.145, 0.035, 0.145, 0.085);
                    auto f3 = std::make_unique<TH1D>(
                        ("h_rob_pull_" + cfg.output_token).c_str(),
                        ";#theta_{p} (deg);Regional-only residual pull",
                        100, pcfg->edges.front(), pcfg->edges.back());
                    f3->SetDirectory(nullptr);
                    f3->SetMinimum(-production_pull_limit());
                    f3->SetMaximum(production_pull_limit());
                    style_production_frame(f3.get(), 0.050, 0.050, 0.041);
                    f3->Draw();
                    rob_owned.emplace_back(std::move(f3));
                    auto z3 = std::make_unique<TLine>(
                        pcfg->edges.front(), 0.0, pcfg->edges.back(), 0.0);
                    z3->SetLineStyle(2); z3->SetLineWidth(2);
                    z3->SetLineColor(kGray + 1); z3->Draw("SAME");
                    rob_owned.emplace_back(std::move(z3));

                    std::array<double, 7> sw{}, swy{};
                    for (const RelativeSlopePoint& pt : full_pts) {
                        if (pt.region_index < 0 || pt.region_index >= 7 ||
                            !(pt.slope_err > 0.0)) continue;
                        const double w = 1.0 / (pt.slope_err * pt.slope_err);
                        sw[pt.region_index] += w;
                        swy[pt.region_index] += w * pt.slope;
                    }
                    std::array<double, 7> mean{};
                    for (int ir = 0; ir < 7; ++ir)
                        mean[ir] = sw[ir] > 0.0 ? swy[ir] / sw[ir] : 0.0;

                    std::array<std::unique_ptr<TGraph>, 7> pg;
                    for (int ir = 0; ir < 7; ++ir) {
                        pg[ir] = std::make_unique<TGraph>();
                        pg[ir]->SetMarkerStyle(photon_region_marker(ir));
                        pg[ir]->SetMarkerColor(photon_region_color(ir));
                        pg[ir]->SetLineColor(photon_region_color(ir));
                        pg[ir]->SetMarkerSize(1.15);
                    }
                    std::array<int,7> np{};
                    for (const RelativeSlopePoint& pt : full_pts) {
                        if (pt.region_index < 0 || pt.region_index >= 7 ||
                            !(pt.slope_err > 0.0)) continue;
                        const double pull =
                            (pt.slope - mean[pt.region_index]) / pt.slope_err;
                        pg[pt.region_index]->SetPoint(
                            np[pt.region_index]++, pt.x, pull);
                    }
                    for (int ir = 0; ir < 7; ++ir) {
                        pg[ir]->Draw("P SAME");
                        rob_owned.emplace_back(std::move(pg[ir]));
                    }
                    auto p3lab = std::make_unique<TLatex>();
                    p3lab->SetNDC(); p3lab->SetTextFont(62); p3lab->SetTextSize(0.047);
                    p3lab->DrawLatex(0.17, 0.88, "Why the full fit looks significant");
                    rob_owned.emplace_back(std::move(p3lab));

                    // Panel 4: full vs omit-first-bin DeltaBIC.
                    crob.cd(4);
                    configure_production_pad((TPad*)gPad, 0.16, 0.035, 0.18, 0.085);
                    const double dbic_full =
                        full_cmp.valid ? full_cmp.bic_m1 - full_cmp.bic_m2 : 0.0;
                    const double dbic_trim =
                        trim_cmp.valid ? trim_cmp.bic_m1 - trim_cmp.bic_m2 : 0.0;
                    const double ymax =
                        std::max(12.0, 1.18 * std::max(dbic_full, dbic_trim));
                    auto f4 = std::make_unique<TH1D>(
                        ("h_rob_bic_" + cfg.output_token).c_str(),
                        ";Fit sample;#DeltaBIC (regional - regional+#theta_{p})",
                        2, 0.5, 2.5);
                    f4->SetDirectory(nullptr);
                    f4->GetXaxis()->SetBinLabel(1, "all #theta_{p} bins");
                    f4->GetXaxis()->SetBinLabel(2, "omit 10-18^{#circ}");
                    f4->SetMinimum(std::min(-5.0, 1.15 * std::min(dbic_full, dbic_trim)));
                    f4->SetMaximum(ymax);
                    style_production_frame(f4.get(), 0.050, 0.050, 0.037);
                    f4->Draw();
                    rob_owned.emplace_back(std::move(f4));
                    auto gb = std::make_unique<TGraph>();
                    gb->SetPoint(0, 1.0, dbic_full);
                    gb->SetPoint(1, 2.0, dbic_trim);
                    gb->SetMarkerStyle(20); gb->SetMarkerSize(1.5);
                    gb->Draw("P SAME");
                    rob_owned.emplace_back(std::move(gb));
                    auto z4 = std::make_unique<TLine>(0.5, 0.0, 2.5, 0.0);
                    z4->SetLineStyle(2); z4->SetLineColor(kGray + 1);
                    z4->Draw("SAME"); rob_owned.emplace_back(std::move(z4));
                    auto info = std::make_unique<TLatex>();
                    info->SetNDC(); info->SetTextFont(42); info->SetTextSize(0.035);
                    std::ostringstream ss1, ss2;
                    ss1 << "full gradient significance = "
                        << std::fixed << std::setprecision(2)
                        << full_cmp.gradient_significance << "#sigma";
                    ss2 << "without first bin = "
                        << std::fixed << std::setprecision(2)
                        << trim_cmp.gradient_significance << "#sigma";
                    info->DrawLatex(0.18, 0.86, ss1.str().c_str());
                    info->DrawLatex(0.18, 0.80, ss2.str().c_str());
                    rob_owned.emplace_back(std::move(info));

                    crob.cd(0);
                    auto title = std::make_unique<TLatex>();
                    title->SetNDC(); title->SetTextAlign(22); title->SetTextFont(62);
                    title->SetTextSize(0.030);
                    title->DrawLatex(
                        0.5, 0.972,
                        "Fa18 Out: robustness of the apparent proton-angle current response");
                    rob_owned.emplace_back(std::move(title));

                    crob.Modified();
                    crob.Update();
                    const std::string path =
                        odir + "/fa18_out_theta_p_lowest_bin_robustness.png";
                    crob.SaveAs(path.c_str());
                    gSystem->CopyFile(
                        path.c_str(),
                        (note_dir + "/" + cfg.output_token +
                         "_fa18_out_theta_p_lowest_bin_robustness.png").c_str(),
                        true);
                }
            }
        }
    }


    // Rank candidate kinematic variables by the improvement of M2 over the
    // region-only model M1.  Positive DeltaBIC means the additional linear
    // dependence is preferred after the parameter penalty.
    std::ofstream ranking_csv(odir + "/regional_kinematic_model_ranking.csv");
    ranking_csv
        << "sample,period,rank,variable,delta_bic_m1_to_m2,"
        << "delta_aic_m1_to_m2,delta_chi2_m1_to_m2,"
        << "gradient,gradient_stat,gradient_significance,m2_chi2,m2_ndf\n";

    for (const std::string& sample : {std::string("data"), std::string("mc")}) {
        if (sample == "mc" && !process_mc) continue;
        for (const std::string& period : periods) {
            std::vector<CandidateModelRow> rows;
            for (const CandidateModelRow& row : model_rows) {
                if (row.sample == sample && row.period == period && row.cmp.valid)
                    rows.push_back(row);
            }
            std::sort(rows.begin(), rows.end(),
                      [](const CandidateModelRow& a, const CandidateModelRow& b) {
                          return (a.cmp.bic_m1 - a.cmp.bic_m2) >
                                 (b.cmp.bic_m1 - b.cmp.bic_m2);
                      });
            for (size_t irank = 0; irank < rows.size(); ++irank) {
                const auto& row = rows[irank];
                ranking_csv << sample << "," << period << "," << (irank + 1) << ","
                            << row.variable << ","
                            << (row.cmp.bic_m1 - row.cmp.bic_m2) << ","
                            << (row.cmp.aic_m1 - row.cmp.aic_m2) << ","
                            << (row.cmp.chi2_m1 - row.cmp.chi2_m2) << ","
                            << row.cmp.gradient << "," << row.cmp.gradient_err << ","
                            << row.cmp.gradient_significance << ","
                            << row.cmp.chi2_m2 << "," << row.cmp.ndf_m2 << "\n";
            }

            if (!rows.empty()) {
                std::cout << "[current_dependence] model selection " << cfg.csv_channel
                          << " " << sample << " " << period
                          << ": best residual variable=" << rows.front().variable
                          << ", DeltaBIC(M1-M2)="
                          << (rows.front().cmp.bic_m1 - rows.front().cmp.bic_m2)
                          << ", gradient significance="
                          << rows.front().cmp.gradient_significance << " sigma"
                          << std::endl;
            }
        }
    }

    // Leave-one-point-out influence study for the three polar-angle
    // candidates.  No event-count threshold is imposed here.  Instead, this
    // reports directly whether a model preference is being driven by one or two
    // individual cells (as occurred in the low-theta_p Fa18 Out edge bins).
    std::ofstream influence_csv(
        odir + "/polar_angle_model_influence.csv");
    influence_csv
        << "sample,period,variable,region,x,total_counts,n_current_points,"
        << "full_delta_bic,loo_delta_bic,delta_bic_influence\n";

    for (const std::string& sample : {std::string("data"), std::string("mc")}) {
        if (sample == "mc" && !process_mc) continue;
        for (const std::string& period : periods) {
            for (const std::string& variable :
                 {std::string("e_theta"), std::string("p_theta"), std::string("g_theta")}) {

                const KinematicVarConfig* vcfg = find_var_config(variable);
                if (!vcfg) continue;

                std::vector<RelativeSlopePoint> pts;
                if (sample == "data") {
                    std::map<std::string, std::vector<DataAgg>> dmap;
                    for (const PhotonRegionSpec& spec : PHOTON_REGION_ORDER)
                        dmap[spec.key] = merged[period][spec.key][variable];
                    pts = relative_slope_points_from_data_theta_bins(dmap, *vcfg);
                } else {
                    std::map<std::string, std::vector<McKinematicBinAgg>> mmap;
                    for (const PhotonRegionSpec& spec : PHOTON_REGION_ORDER)
                        mmap[spec.key] = merged_mc[period][spec.key][variable];
                    pts = relative_slope_points_from_mc_theta_bins(mmap, *vcfg);
                }

                const RelativeSlopeModelComparison full =
                    compare_relative_slope_models(pts);
                if (!full.valid) continue;
                const double full_dbic = full.bic_m1 - full.bic_m2;

                for (size_t i = 0; i < pts.size(); ++i) {
                    std::vector<RelativeSlopePoint> loo;
                    loo.reserve(pts.size() - 1);
                    for (size_t j = 0; j < pts.size(); ++j)
                        if (j != i) loo.push_back(pts[j]);

                    const RelativeSlopeModelComparison cmp =
                        compare_relative_slope_models(loo);
                    if (!cmp.valid) continue;
                    const double loo_dbic = cmp.bic_m1 - cmp.bic_m2;
                    const RelativeSlopePoint& pt = pts[i];
                    const std::string region =
                        (pt.region_index >= 0 && pt.region_index < 7)
                            ? PHOTON_REGION_ORDER[pt.region_index].key
                            : "unknown";

                    influence_csv
                        << sample << "," << period << "," << variable << ","
                        << region << "," << pt.x << ","
                        << pt.total_counts << "," << pt.n_current_points << ","
                        << full_dbic << "," << loo_dbic << ","
                        << (full_dbic - loo_dbic) << "\n";
                }
            }
        }
    }

    std::ofstream polar_rank_csv(
        odir + "/polar_angle_model_ranking.csv");
    polar_rank_csv
        << "sample,period,rank,variable,delta_bic_m1_to_m2,delta_aic_m1_to_m2,"
        << "delta_chi2_m1_to_m2,gradient,gradient_stat,gradient_significance,"
        << "m2_chi2,m2_ndf\n";

    for (const std::string& sample : {std::string("data"), std::string("mc")}) {
        if (sample == "mc" && !process_mc) continue;
        for (const std::string& period : periods) {
            std::vector<const CandidateModelRow*> candidates;
            for (const CandidateModelRow& row : model_rows) {
                if (row.sample == sample && row.period == period &&
                    row.cmp.valid && is_polar_angle_model_variable(row.variable)) {
                    candidates.push_back(&row);
                }
            }
            std::sort(candidates.begin(), candidates.end(),
                      [](const CandidateModelRow* a, const CandidateModelRow* b) {
                          return (a->cmp.bic_m1 - a->cmp.bic_m2) >
                                 (b->cmp.bic_m1 - b->cmp.bic_m2);
                      });
            for (size_t i = 0; i < candidates.size(); ++i) {
                const CandidateModelRow& row = *candidates[i];
                polar_rank_csv
                    << sample << "," << period << "," << (i + 1) << ","
                    << row.variable << ","
                    << row.cmp.bic_m1 - row.cmp.bic_m2 << ","
                    << row.cmp.aic_m1 - row.cmp.aic_m2 << ","
                    << row.cmp.chi2_m1 - row.cmp.chi2_m2 << ","
                    << row.cmp.gradient << "," << row.cmp.gradient_err << ","
                    << row.cmp.gradient_significance << ","
                    << row.cmp.chi2_m2 << "," << row.cmp.ndf_m2 << "\n";
            }
        }
    }

    // Dedicated Fa18-Out theta_p edge-bin stress test.  This is not a
    // statistical-count cut: the diagnostic simply repeats the identical
    // nested M1/M2 comparison after removing the configured lowest theta_p bin
    // (10--18 deg) from every photon region.  The purpose is to expose whether
    // the apparent theta_p preference is a stable trend across the populated
    // acceptance or is driven by that kinematic edge.
    if (cfg.channel == Channel::DVCS) {
        const KinematicVarConfig* ptheta_cfg = find_var_config("p_theta");
        if (ptheta_cfg && ptheta_cfg->edges.size() >= 2) {
            std::map<std::string, std::vector<DataAgg>> dmap;
            for (const PhotonRegionSpec& spec : PHOTON_REGION_ORDER)
                dmap[spec.key] = merged["Fa18 Out"][spec.key]["p_theta"];

            const std::vector<RelativeSlopePoint> full_pts =
                relative_slope_points_from_data_theta_bins(dmap, *ptheta_cfg);
            std::vector<RelativeSlopePoint> trimmed_pts;
            const double first_bin_high = ptheta_cfg->edges[1];
            for (const RelativeSlopePoint& pt : full_pts) {
                if (pt.x >= first_bin_high) trimmed_pts.push_back(pt);
            }

            const RelativeSlopeModelComparison full_cmp =
                compare_relative_slope_models(full_pts);
            const RelativeSlopeModelComparison trim_cmp =
                compare_relative_slope_models(trimmed_pts);
            const PooledRegionThetaSlopeFit full_fit =
                fit_pooled_region_theta_relative_slopes(full_pts);
            const PooledRegionThetaSlopeFit trim_fit =
                fit_pooled_region_theta_relative_slopes(trimmed_pts);

            std::ofstream stress_csv(
                odir + "/fa18_out_theta_p_low_edge_stress_test.csv");
            stress_csv
                << "selection,npoints,m1_chi2,m1_ndf,m1_chi2_per_ndf,"
                << "m2_chi2,m2_ndf,m2_chi2_per_ndf,delta_bic_m1_to_m2,"
                << "gradient_per_nA_per_deg,gradient_stat,gradient_significance\n";

            auto write_stress_row = [&](const std::string& label,
                                        const RelativeSlopeModelComparison& cmp,
                                        const PooledRegionThetaSlopeFit& fit) {
                stress_csv << label << "," << cmp.npoints << ","
                           << cmp.chi2_m1 << "," << cmp.ndf_m1 << ","
                           << (cmp.ndf_m1 > 0 ? cmp.chi2_m1 / cmp.ndf_m1
                                             : std::numeric_limits<double>::quiet_NaN()) << ","
                           << cmp.chi2_m2 << "," << cmp.ndf_m2 << ","
                           << (cmp.ndf_m2 > 0 ? cmp.chi2_m2 / cmp.ndf_m2
                                             : std::numeric_limits<double>::quiet_NaN()) << ","
                           << (cmp.bic_m1 - cmp.bic_m2) << ","
                           << fit.theta_gradient << "," << fit.theta_gradient_err << ","
                           << ((fit.valid && fit.theta_gradient_err > 0.0)
                                   ? std::fabs(fit.theta_gradient) / fit.theta_gradient_err
                                   : std::numeric_limits<double>::quiet_NaN()) << "\n";
            };
            if (full_cmp.valid && full_fit.valid)
                write_stress_row("all_theta_p_bins", full_cmp, full_fit);
            if (trim_cmp.valid && trim_fit.valid)
                write_stress_row("exclude_10_to_18_deg", trim_cmp, trim_fit);

            std::ofstream point_csv(
                odir + "/fa18_out_theta_p_low_edge_points.csv");
            point_csv << "region,theta_p_deg,relative_slope_per_nA,relative_slope_stat,"
                      << "total_counts,n_current_points,is_lowest_theta_bin\n";
            for (const RelativeSlopePoint& pt : full_pts) {
                const std::string region =
                    (pt.region_index >= 0 && pt.region_index < 7)
                        ? PHOTON_REGION_ORDER[pt.region_index].key : "unknown";
                point_csv << region << "," << pt.x << "," << pt.slope << ","
                          << pt.slope_err << "," << pt.total_counts << ","
                          << pt.n_current_points << ","
                          << (pt.x < first_bin_high ? 1 : 0) << "\n";
            }

            TCanvas cstress(
                ("c_fa18_out_theta_p_edge_stress_" + cfg.output_token).c_str(),
                "", 2250, 720);
            cstress.Divide(3, 1, 0.002, 0.002);
            std::vector<std::unique_ptr<TObject>> owned;

            // Panel 1: measured current-response slopes.  The lowest theta_p
            // bin is intentionally drawn open so its leverage is obvious.
            cstress.cd(1);
            configure_production_pad((TPad*)gPad, 0.15, 0.035, 0.15, 0.10);
            double central_lo = std::numeric_limits<double>::infinity();
            double central_hi = -std::numeric_limits<double>::infinity();
            for (const RelativeSlopePoint& pt : trimmed_pts) {
                if (!std::isfinite(pt.slope)) continue;
                central_lo = std::min(central_lo, 100.0 * pt.slope);
                central_hi = std::max(central_hi, 100.0 * pt.slope);
            }
            if (!std::isfinite(central_lo) || !std::isfinite(central_hi) ||
                !(central_hi > central_lo)) {
                central_lo = -1.0;
                central_hi = 0.5;
            }
            double slope_span = std::max(0.25, central_hi - central_lo);
            const double slope_ymin = central_lo - 0.25 * slope_span;
            const double slope_ymax = central_hi + 0.25 * slope_span;
            auto slope_frame = std::make_unique<TH1D>(
                ("h_fa18_out_ptheta_stress_points_" + cfg.output_token).c_str(),
                ";#theta_{p} (deg);Relative current slope (%/nA)",
                100, ptheta_cfg->edges.front(), ptheta_cfg->edges.back());
            slope_frame->SetDirectory(nullptr);
            slope_frame->SetMinimum(slope_ymin);
            slope_frame->SetMaximum(slope_ymax);
            style_production_frame(slope_frame.get(), 0.052, 0.052, 0.041);
            slope_frame->Draw();
            owned.emplace_back(std::move(slope_frame));

            std::array<std::unique_ptr<TGraphErrors>, 7> good_graphs;
            std::array<std::unique_ptr<TGraphErrors>, 7> edge_graphs;
            std::array<int, 7> ngood{};
            std::array<int, 7> nedge{};
            for (int ir = 0; ir < 7; ++ir) {
                good_graphs[ir] = std::make_unique<TGraphErrors>();
                edge_graphs[ir] = std::make_unique<TGraphErrors>();
                good_graphs[ir]->SetMarkerStyle(photon_region_marker(ir));
                edge_graphs[ir]->SetMarkerStyle(24 + (ir % 4));
                good_graphs[ir]->SetMarkerColor(photon_region_color(ir));
                edge_graphs[ir]->SetMarkerColor(photon_region_color(ir));
                good_graphs[ir]->SetLineColor(photon_region_color(ir));
                edge_graphs[ir]->SetLineColor(photon_region_color(ir));
                good_graphs[ir]->SetMarkerSize(1.15);
                edge_graphs[ir]->SetMarkerSize(1.15);
            }
            for (const RelativeSlopePoint& pt : full_pts) {
                if (pt.region_index < 0 || pt.region_index >= 7) continue;
                TGraphErrors* g = pt.x < first_bin_high
                    ? edge_graphs[pt.region_index].get()
                    : good_graphs[pt.region_index].get();
                int& n = pt.x < first_bin_high
                    ? nedge[pt.region_index] : ngood[pt.region_index];
                g->SetPoint(n, pt.x, 100.0 * pt.slope);
                g->SetPointError(n, 0.0, 100.0 * pt.slope_err);
                ++n;
            }
            for (int ir = 0; ir < 7; ++ir) {
                good_graphs[ir]->Draw("P SAME");
                edge_graphs[ir]->Draw("P SAME");
                owned.emplace_back(std::move(good_graphs[ir]));
                owned.emplace_back(std::move(edge_graphs[ir]));
            }
            auto edge_line = std::make_unique<TLine>(
                first_bin_high, slope_ymin, first_bin_high, slope_ymax);
            edge_line->SetLineStyle(2);
            edge_line->SetLineWidth(2);
            edge_line->SetLineColor(kGray + 1);
            edge_line->Draw("SAME");
            owned.emplace_back(std::move(edge_line));
            auto edge_note = std::make_unique<TLatex>();
            edge_note->SetNDC();
            edge_note->SetTextFont(42);
            edge_note->SetTextSize(0.035);
            edge_note->DrawLatex(0.18, 0.84, "open: 10^{#circ}<#theta_{p}<18^{#circ} edge bin");
            owned.emplace_back(std::move(edge_note));

            // Panel 2: model preference before/after removing only that edge bin.
            cstress.cd(2);
            configure_production_pad((TPad*)gPad, 0.16, 0.035, 0.18, 0.10);
            auto bic_frame = std::make_unique<TH1D>(
                ("h_fa18_out_ptheta_stress_bic_" + cfg.output_token).c_str(),
                ";Fit sample;#DeltaBIC = BIC_{regional}-BIC_{regional+#theta_{p}}",
                2, 0.5, 2.5);
            bic_frame->SetDirectory(nullptr);
            bic_frame->GetXaxis()->SetBinLabel(1, "all bins");
            bic_frame->GetXaxis()->SetBinLabel(2, "exclude 10-18^{#circ}");
            const double full_dbic = full_cmp.valid ? full_cmp.bic_m1 - full_cmp.bic_m2 : 0.0;
            const double trim_dbic = trim_cmp.valid ? trim_cmp.bic_m1 - trim_cmp.bic_m2 : 0.0;
            const double bic_max = std::max(12.0, 1.15 * std::max(full_dbic, trim_dbic));
            bic_frame->SetMinimum(std::min(-5.0, 1.15 * std::min(full_dbic, trim_dbic)));
            bic_frame->SetMaximum(bic_max);
            style_production_frame(bic_frame.get(), 0.052, 0.052, 0.041);
            bic_frame->Draw();
            owned.emplace_back(std::move(bic_frame));
            auto bic_zero = std::make_unique<TLine>(0.5, 0.0, 2.5, 0.0);
            bic_zero->SetLineStyle(2);
            bic_zero->SetLineWidth(2);
            bic_zero->SetLineColor(kGray + 1);
            bic_zero->Draw("SAME");
            owned.emplace_back(std::move(bic_zero));
            auto bic_graph = std::make_unique<TGraph>();
            bic_graph->SetPoint(0, 1.0, full_dbic);
            bic_graph->SetPoint(1, 2.0, trim_dbic);
            bic_graph->SetMarkerStyle(20);
            bic_graph->SetMarkerSize(1.5);
            bic_graph->Draw("P SAME");
            owned.emplace_back(std::move(bic_graph));

            // Panel 3: absolute closure quality.  If the regional-only result
            // becomes satisfactory once the edge bin is removed, there is no
            // motivation to promote theta_p into the production correction.
            cstress.cd(3);
            configure_production_pad((TPad*)gPad, 0.16, 0.035, 0.18, 0.10);
            auto chi_frame = std::make_unique<TH1D>(
                ("h_fa18_out_ptheta_stress_chi2_" + cfg.output_token).c_str(),
                ";Fit sample;#chi^{2}/ndf",
                4, 0.5, 4.5);
            chi_frame->SetDirectory(nullptr);
            chi_frame->GetXaxis()->SetBinLabel(1, "M1 all");
            chi_frame->GetXaxis()->SetBinLabel(2, "M2 all");
            chi_frame->GetXaxis()->SetBinLabel(3, "M1 no edge");
            chi_frame->GetXaxis()->SetBinLabel(4, "M2 no edge");
            std::array<double, 4> chi_vals = {
                full_cmp.valid && full_cmp.ndf_m1 > 0 ? full_cmp.chi2_m1 / full_cmp.ndf_m1 : 0.0,
                full_cmp.valid && full_cmp.ndf_m2 > 0 ? full_cmp.chi2_m2 / full_cmp.ndf_m2 : 0.0,
                trim_cmp.valid && trim_cmp.ndf_m1 > 0 ? trim_cmp.chi2_m1 / trim_cmp.ndf_m1 : 0.0,
                trim_cmp.valid && trim_cmp.ndf_m2 > 0 ? trim_cmp.chi2_m2 / trim_cmp.ndf_m2 : 0.0
            };
            double chi_max = 2.0;
            for (double x : chi_vals) chi_max = std::max(chi_max, 1.15 * x);
            chi_frame->SetMinimum(0.0);
            chi_frame->SetMaximum(chi_max);
            style_production_frame(chi_frame.get(), 0.052, 0.052, 0.041);
            chi_frame->Draw();
            owned.emplace_back(std::move(chi_frame));
            auto chi_one = std::make_unique<TLine>(0.5, 1.0, 4.5, 1.0);
            chi_one->SetLineStyle(2);
            chi_one->SetLineWidth(2);
            chi_one->SetLineColor(kGray + 1);
            chi_one->Draw("SAME");
            owned.emplace_back(std::move(chi_one));
            auto chi_graph = std::make_unique<TGraph>();
            for (int i = 0; i < 4; ++i) chi_graph->SetPoint(i, i + 1.0, chi_vals[i]);
            chi_graph->SetMarkerStyle(20);
            chi_graph->SetMarkerSize(1.4);
            chi_graph->Draw("P SAME");
            owned.emplace_back(std::move(chi_graph));

            cstress.cd(0);
            auto stress_title = std::make_unique<TLatex>();
            stress_title->SetNDC();
            stress_title->SetTextAlign(22);
            stress_title->SetTextFont(62);
            stress_title->SetTextSize(0.031);
            stress_title->DrawLatex(
                0.5, 0.968,
                "Fa18 Out: #theta_{p} response-model edge-bin stress test");
            owned.emplace_back(std::move(stress_title));
            cstress.Modified();
            cstress.Update();

            const std::string stress_path =
                odir + "/fa18_out_theta_p_low_edge_stress_test.png";
            cstress.SaveAs(stress_path.c_str());
            gSystem->CopyFile(
                stress_path.c_str(),
                (note_dir + "/" + cfg.output_token +
                 "_fa18_out_theta_p_low_edge_stress_test.png").c_str(),
                true);

            std::cout << "[current_dependence] Fa18 Out theta_p edge-bin stress test: "
                      << "DeltaBIC full=" << full_dbic
                      << ", without 10-18 deg=" << trim_dbic
                      << "; M1 chi2/ndf full="
                      << (full_cmp.ndf_m1 > 0 ? full_cmp.chi2_m1/full_cmp.ndf_m1 : 0.0)
                      << ", without edge="
                      << (trim_cmp.ndf_m1 > 0 ? trim_cmp.chi2_m1/trim_cmp.ndf_m1 : 0.0)
                      << std::endl;
        }
    }

    // Polar-angle-only model-selection canvas used for response-model
    // decisions.  The broader seven-variable plot remains available as a
    // kinematic correlation diagnostic, but only theta_e, theta_p and
    // theta_gamma are treated as candidate detector-response coordinates.
    TCanvas cpolar(("c_polar_angle_model_selection_" + cfg.output_token).c_str(),
                   "", 1800, 1200);
    cpolar.Divide(2, 2, 0.002, 0.002);
    std::vector<std::unique_ptr<TObject>> polar_owned;
    const std::vector<std::string> polar_vars = {
        "e_theta", "p_theta", "g_theta"
    };

    for (size_t ip = 0; ip < periods.size(); ++ip) {
        cpolar.cd((int)ip + 1);
        configure_production_pad((TPad*)gPad, 0.145, 0.035, 0.16, 0.075);

        auto frame = std::make_unique<TH1D>(
            ("h_polar_model_selection_" + cfg.output_token + "_" +
             std::to_string(ip)).c_str(),
            ";Candidate polar angle;#DeltaBIC = BIC_{regional}-BIC_{regional+angle}",
            3, 0.5, 3.5);
        frame->SetDirectory(nullptr);
        for (int iv = 0; iv < 3; ++iv)
            frame->GetXaxis()->SetBinLabel(
                iv + 1, pretty_kinematic_key(polar_vars[iv]).c_str());

        double ymin = -5.0;
        double ymax = 5.0;
        std::array<double, 3> values = {
            std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN()
        };
        for (int iv = 0; iv < 3; ++iv) {
            for (const CandidateModelRow& row : model_rows) {
                if (row.sample == "data" && row.period == periods[ip] &&
                    row.variable == polar_vars[iv] && row.cmp.valid) {
                    values[iv] = row.cmp.bic_m1 - row.cmp.bic_m2;
                    ymin = std::min(ymin, values[iv]);
                    ymax = std::max(ymax, values[iv]);
                }
            }
        }
        const double span = std::max(5.0, ymax - ymin);
        frame->SetMinimum(std::min(-5.0, ymin - 0.12 * span));
        frame->SetMaximum(std::max(5.0, ymax + 0.15 * span));
        style_production_frame(frame.get(), 0.050, 0.050, 0.040);
        frame->Draw();
        polar_owned.emplace_back(std::move(frame));

        auto zero = std::make_unique<TLine>(0.5, 0.0, 3.5, 0.0);
        zero->SetLineStyle(2);
        zero->SetLineWidth(2);
        zero->SetLineColor(kGray + 1);
        zero->Draw("SAME");
        polar_owned.emplace_back(std::move(zero));

        auto graph = std::make_unique<TGraph>();
        int ng = 0;
        for (int iv = 0; iv < 3; ++iv)
            if (std::isfinite(values[iv]))
                graph->SetPoint(ng++, iv + 1.0, values[iv]);
        graph->SetMarkerStyle(20);
        graph->SetMarkerSize(1.4);
        graph->SetLineWidth(2);
        graph->Draw("P SAME");
        polar_owned.emplace_back(std::move(graph));

        polar_owned.emplace_back(
            make_period_label(periods[ip], 0.94, 0.86, 33, 0.052));
    }

    cpolar.cd(0);
    auto ptitle = std::make_unique<TLatex>();
    ptitle->SetNDC();
    ptitle->SetTextAlign(22);
    ptitle->SetTextFont(62);
    ptitle->SetTextSize(0.029);
    ptitle->DrawLatex(
        0.5, 0.975,
        (cfg.title + ": polar-angle current-response model selection").c_str());
    polar_owned.emplace_back(std::move(ptitle));

    cpolar.Modified();
    cpolar.Update();
    const std::string polar_path =
        odir + "/polar_angle_model_selection_delta_bic.png";
    cpolar.SaveAs(polar_path.c_str());
    gSystem->CopyFile(
        polar_path.c_str(),
        (note_dir + "/" + cfg.output_token +
         "_polar_angle_model_selection_delta_bic.png").c_str(),
        true);

    // Production-quality four-panel DATA ranking plot.  Each panel is
    // explicitly labeled with its run period.  Larger positive DeltaBIC favors
    // regional + variable over the region-only model.
    TCanvas csel(("c_region_model_selection_" + cfg.output_token).c_str(), "", 1800, 1200);
    csel.Divide(2, 2, 0.002, 0.002);
    std::vector<std::unique_ptr<TObject>> selection_owned;

    for (size_t ip = 0; ip < periods.size(); ++ip) {
        csel.cd((int)ip + 1);
        configure_production_pad((TPad*)gPad, 0.145, 0.035, 0.18, 0.075);

        std::vector<CandidateModelRow> rows;
        for (const CandidateModelRow& row : model_rows) {
            if (row.sample == "data" && row.period == periods[ip] && row.cmp.valid)
                rows.push_back(row);
        }

        auto frame = std::make_unique<TH1D>(
            ("h_model_selection_" + cfg.output_token + "_" + std::to_string(ip)).c_str(),
            ";Candidate residual variable;#DeltaBIC = BIC_{regional}-BIC_{regional+var}",
            (int)vars.size(), 0.5, vars.size() + 0.5);
        frame->SetDirectory(nullptr);
        for (size_t iv = 0; iv < vars.size(); ++iv)
            frame->GetXaxis()->SetBinLabel(
                (int)iv + 1, pretty_kinematic_key(vars[iv].key).c_str());

        double ymin = 0.0;
        double ymax = 0.0;
        bool have = false;
        for (const CandidateModelRow& row : rows) {
            const double d = row.cmp.bic_m1 - row.cmp.bic_m2;
            if (!std::isfinite(d)) continue;
            ymin = have ? std::min(ymin, d) : d;
            ymax = have ? std::max(ymax, d) : d;
            have = true;
        }
        if (!have) { ymin = -5.0; ymax = 5.0; }
        const double span = std::max(5.0, ymax - ymin);
        frame->SetMinimum(std::min(-5.0, ymin - 0.12 * span));
        frame->SetMaximum(std::max(5.0, ymax + 0.15 * span));
        style_production_frame(frame.get(), 0.050, 0.050, 0.041);
        frame->GetXaxis()->LabelsOption("h");
        frame->Draw();
        selection_owned.emplace_back(std::move(frame));

        auto zero = std::make_unique<TLine>(0.5, 0.0, vars.size() + 0.5, 0.0);
        zero->SetLineStyle(2);
        zero->SetLineWidth(2);
        zero->SetLineColor(kGray + 1);
        zero->Draw("SAME");
        selection_owned.emplace_back(std::move(zero));

        // Conventional BIC evidence guides: +6 strong, +10 very strong.
        for (double threshold : {6.0, 10.0}) {
            auto guide = std::make_unique<TLine>(
                0.5, threshold, vars.size() + 0.5, threshold);
            guide->SetLineStyle(threshold == 6.0 ? 3 : 7);
            guide->SetLineWidth(1);
            guide->SetLineColor(kGray + 1);
            guide->Draw("SAME");
            selection_owned.emplace_back(std::move(guide));
        }

        auto graph = std::make_unique<TGraph>();
        int ng = 0;
        for (size_t iv = 0; iv < vars.size(); ++iv) {
            for (const CandidateModelRow& row : rows) {
                if (row.variable != vars[iv].key) continue;
                const double d = row.cmp.bic_m1 - row.cmp.bic_m2;
                if (!std::isfinite(d)) continue;
                graph->SetPoint(ng++, (double)iv + 1.0, d);
            }
        }
        graph->SetMarkerStyle(20);
        graph->SetMarkerSize(1.35);
        graph->SetLineWidth(2);
        graph->Draw("P SAME");
        selection_owned.emplace_back(std::move(graph));

        selection_owned.emplace_back(make_period_label(periods[ip], 0.94, 0.86, 33, 0.052));
    }

    csel.cd(0);
    auto stitle = std::make_unique<TLatex>();
    stitle->SetNDC();
    stitle->SetTextAlign(22);
    stitle->SetTextFont(62);
    stitle->SetTextSize(0.030);
    stitle->DrawLatex(
        0.5, 0.972,
        (cfg.title + ": residual current-response model selection after photon-region correction").c_str());
    selection_owned.emplace_back(std::move(stitle));
    csel.Modified();
    csel.Update();

    const std::string selection_path =
        odir + "/regional_kinematic_model_selection_delta_bic.png";
    csel.SaveAs(selection_path.c_str());
    gSystem->CopyFile(
        selection_path.c_str(),
        (note_dir + "/" + cfg.output_token + "_regional_kinematic_model_selection_delta_bic.png").c_str(),
        true);

    std::cout << "[current_dependence] Region-conditioned kinematic diagnostics written under "
              << odir << std::endl;
    return production_e_theta_models;
}

static void write_photon_region_summary_csv(const std::string& path,
                                            const PhotonRegionResults& results,
                                            const std::vector<PeriodResult>& integrated_results) {
    const size_t slash_pos = path.find_last_of('/');
    if (slash_pos != std::string::npos) {
        mkdir_p(path.substr(0, slash_pos));
    }

    std::ofstream fout(path);
    if (!fout.is_open()) {
        fatal("[current_dependence] FATAL: cannot write photon-region diagnostic CSV: " + path);
    }

    fout << "period,photon_region,data_factor,data_factor_stat,mc_factor,mc_factor_stat,"
         << "data_over_mc,data_over_mc_stat,"
         << "data_over_mc_relative_to_integrated,data_over_mc_relative_to_integrated_stat,"
         << "data_slope_percent_per_nA,data_slope_percent_per_nA_stat,"
         << "mc_slope_percent_per_nA,mc_slope_percent_per_nA_stat,"
         << "integrated_data_factor,integrated_mc_factor,integrated_data_over_mc\n";

    for (const std::string& period : PERIOD_ORDER) {
        const PeriodResult* integrated = find_period_result(integrated_results, period);

        for (const PhotonRegionSpec& spec : PHOTON_REGION_ORDER) {
            auto it_region = results.find(spec.key);
            if (it_region == results.end()) {
                continue;
            }
            const PeriodResult* r = find_period_result(it_region->second, period);
            if (!r) {
                continue;
            }

            double data_over_mc = std::numeric_limits<double>::quiet_NaN();
            double data_over_mc_err = std::numeric_limits<double>::quiet_NaN();
            double int_data_over_mc = std::numeric_limits<double>::quiet_NaN();
            double relative_ratio = std::numeric_limits<double>::quiet_NaN();
            double relative_ratio_err = std::numeric_limits<double>::quiet_NaN();

            if (std::isfinite(r->data_factor) && r->data_factor > 0.0 &&
                std::isfinite(r->mc_factor) && r->mc_factor > 0.0) {
                data_over_mc = r->data_factor / r->mc_factor;
                const double rel2 =
                    std::pow(r->data_factor_err / r->data_factor, 2) +
                    std::pow(r->mc_factor_err / r->mc_factor, 2);
                data_over_mc_err = std::fabs(data_over_mc) *
                                   std::sqrt(std::max(0.0, rel2));
            }

            if (integrated &&
                std::isfinite(integrated->data_factor) && integrated->data_factor > 0.0 &&
                std::isfinite(integrated->mc_factor) && integrated->mc_factor > 0.0) {
                int_data_over_mc = integrated->data_factor / integrated->mc_factor;
                if (std::isfinite(data_over_mc) && int_data_over_mc > 0.0) {
                    relative_ratio = data_over_mc / int_data_over_mc;
                    const double int_rel2 =
                        std::pow(integrated->data_factor_err / integrated->data_factor, 2) +
                        std::pow(integrated->mc_factor_err / integrated->mc_factor, 2);
                    const double reg_rel2 =
                        std::pow(r->data_factor_err / r->data_factor, 2) +
                        std::pow(r->mc_factor_err / r->mc_factor, 2);
                    relative_ratio_err = std::fabs(relative_ratio) *
                        std::sqrt(std::max(0.0, reg_rel2 + int_rel2));
                }
            }

            fout << period << "," << spec.label << ","
                 << std::setprecision(12) << r->data_factor << ","
                 << std::setprecision(12) << r->data_factor_err << ","
                 << std::setprecision(12) << r->mc_factor << ","
                 << std::setprecision(12) << r->mc_factor_err << ","
                 << std::setprecision(12) << data_over_mc << ","
                 << std::setprecision(12) << data_over_mc_err << ","
                 << std::setprecision(12) << relative_ratio << ","
                 << std::setprecision(12) << relative_ratio_err << ","
                 << std::setprecision(12) << fit_percent_slope(r->data_fit) << ","
                 << std::setprecision(12) << fit_percent_slope_err(r->data_fit) << ","
                 << std::setprecision(12) << fit_percent_slope(r->mc_fit) << ","
                 << std::setprecision(12) << fit_percent_slope_err(r->mc_fit) << ","
                 << std::setprecision(12) << (integrated ? integrated->data_factor : std::numeric_limits<double>::quiet_NaN()) << ","
                 << std::setprecision(12) << (integrated ? integrated->mc_factor : std::numeric_limits<double>::quiet_NaN()) << ","
                 << std::setprecision(12) << int_data_over_mc << "\n";
        }
    }
}

static void draw_photon_region_response_canvas(const std::string& out_path,
                                               const PhotonRegionResults& results,
                                               bool data_sample,
                                               bool hide_sp19_inb_from_all_period_plots) {
    mkdir_p(out_path.substr(0, out_path.find_last_of('/')));

    TCanvas c(data_sample ? "c_photon_region_data_response" : "c_photon_region_mc_response",
              data_sample ? "DATA current dependence by photon region" : "MC current dependence by photon region",
              1800, 1200);
    c.Divide(3, 2, 0.002, 0.002);

    const std::vector<std::string> panel_order = {
        "Sp18 Inb", "Sp18 Out", "", "Fa18 Inb", "Fa18 Out", "Sp19 Inb"
    };

    std::vector<TObject*> owned;

    for (int ipad = 0; ipad < 6; ++ipad) {
        c.cd(ipad + 1);
        gPad->SetGridx();
        gPad->SetGridy();

        const std::string period = panel_order[ipad];
        if (period.empty()) {
            TLatex lat;
            lat.SetNDC();
            lat.SetTextSize(0.060);
            lat.DrawLatex(0.10, 0.88, "Photon region");
            for (int ir = 0; ir < (int)PHOTON_REGION_ORDER.size(); ++ir) {
                const double y = 0.76 - 0.085 * ir;
                TLine* line = new TLine(0.14, y, 0.28, y);
                line->SetNDC(true);
                line->SetLineColor(photon_region_color(ir));
                line->SetLineWidth(3);
                line->Draw();
                owned.push_back(line);
                lat.SetTextSize(0.050);
                lat.DrawLatex(0.32, y - 0.018, PHOTON_REGION_ORDER[ir].label.c_str());
            }
            continue;
        }

        TH1D* frame = new TH1D(("h_photon_region_response_" + std::to_string(ipad) + (data_sample ? "_d" : "_m")).c_str(),
                               (period + ";Beam current (nA);Efficiency relative to fitted 0 nA (%)").c_str(),
                               100, 0.0, 80.0);
        frame->SetMinimum(40.0);
        frame->SetMaximum(125.0);
        frame->SetStats(0);
        frame->Draw("AXIS");
        owned.push_back(frame);

        for (int ir = 0; ir < (int)PHOTON_REGION_ORDER.size(); ++ir) {
            auto it = results.find(PHOTON_REGION_ORDER[ir].key);
            if (it == results.end()) {
                continue;
            }
            const PeriodResult* r = find_period_result(it->second, period);
            if (!r) {
                continue;
            }

            const FitResult& fit = data_sample ? r->data_fit : r->mc_fit;
            const std::vector<CurrentPoint>& points = data_sample ? r->data_points : r->mc_points;
            if (!std::isfinite(fit.b) || fit.b == 0.0) {
                continue;
            }

            const int color = photon_region_color(ir);
            TGraphErrors* pts = make_percent_points_graph_style(points, fit, color, photon_region_marker(ir), 0.85);
            TGraph* line = make_percent_fit_line_style(fit, color, 1, 0.0, 80.0);
            if (line) {
                line->SetLineWidth(2);
                line->Draw("L SAME");
                owned.push_back(line);
            }
            if (pts) {
                pts->Draw("PE SAME");
                owned.push_back(pts);
            }
        }

        TLatex* period_label = new TLatex();
        period_label->SetNDC(true);
        period_label->SetTextAlign(33);
        period_label->SetTextFont(62);
        period_label->SetTextSize(0.040);
        period_label->DrawLatex(0.94, 0.88, period.c_str());
        owned.push_back(period_label);

        if (period == "Sp19 Inb" && hide_sp19_inb_from_all_period_plots) {
            TLatex lat;
            lat.SetNDC(true);
            lat.SetTextSize(0.028);
            lat.SetTextColor(kRed + 2);
            lat.DrawLatex(0.13, 0.88, "Direct Sp19 scan shown; nominal DATA response transferred from Fa18 Inb");
        }
    }

    c.SaveAs(out_path.c_str());
    for (TObject* obj : owned) {
        delete obj;
    }
}

static void draw_photon_region_factor_canvas(const std::string& out_path,
                                             const PhotonRegionResults& results,
                                             const std::vector<PeriodResult>& integrated_results,
                                             bool data_sample,
                                             bool hide_sp19_inb_from_all_period_plots) {
    mkdir_p(out_path.substr(0, out_path.find_last_of('/')));

    TCanvas c(data_sample ? "c_photon_region_data_factor" : "c_photon_region_mc_factor",
              data_sample ? "DATA current-efficiency factor by photon region" : "MC current-efficiency factor by photon region",
              1800, 1200);
    c.Divide(3, 2, 0.002, 0.002);

    const std::vector<std::string> panel_order = {
        "Sp18 Inb", "Sp18 Out", "", "Fa18 Inb", "Fa18 Out", "Sp19 Inb"
    };

    std::vector<TObject*> owned;

    for (int ipad = 0; ipad < 6; ++ipad) {
        c.cd(ipad + 1);
        gPad->SetGridy();

        const std::string period = panel_order[ipad];
        if (period.empty()) {
            TLatex lat;
            lat.SetNDC();
            lat.SetTextSize(0.055);
            lat.DrawLatex(0.10, 0.84, data_sample ? "DATA factor" : "MC factor");
            lat.SetTextSize(0.045);
            lat.DrawLatex(0.10, 0.72, "Points: measured FT/S1...S6 current-response factors");
            lat.DrawLatex(0.10, 0.63, "Dashed line: existing integrated factor");
            lat.DrawLatex(0.10, 0.51, "Nominal regional response is applied event-by-event downstream");
            continue;
        }

        TH1D* frame = new TH1D(("h_photon_region_factor_" + std::to_string(ipad) + (data_sample ? "_d" : "_m")).c_str(),
                               (period + ";Photon region;Current-efficiency factor").c_str(),
                               7, -0.5, 6.5);
        frame->SetStats(0);
        frame->SetMinimum(0.45);
        frame->SetMaximum(1.15);
        for (int ir = 0; ir < 7; ++ir) {
            frame->GetXaxis()->SetBinLabel(ir + 1, PHOTON_REGION_ORDER[ir].label.c_str());
        }
        frame->Draw("AXIS");
        owned.push_back(frame);

        TGraphErrors* graph = new TGraphErrors();
        int ip = 0;
        for (int ir = 0; ir < 7; ++ir) {
            auto it = results.find(PHOTON_REGION_ORDER[ir].key);
            if (it == results.end()) continue;
            const PeriodResult* r = find_period_result(it->second, period);
            if (!r) continue;

            const double value = data_sample ? r->data_factor : r->mc_factor;
            const double err = data_sample ? r->data_factor_err : r->mc_factor_err;
            if (!std::isfinite(value)) continue;
            graph->SetPoint(ip, (double)ir, value);
            graph->SetPointError(ip, 0.0, std::isfinite(err) ? err : 0.0);
            ++ip;
        }
        graph->SetMarkerStyle(20);
        graph->SetMarkerSize(1.2);
        graph->SetLineWidth(2);
        graph->Draw("PE SAME");
        owned.push_back(graph);

        const PeriodResult* integrated = find_period_result(integrated_results, period);
        if (integrated) {
            const double baseline = data_sample ? integrated->data_factor : integrated->mc_factor;
            if (std::isfinite(baseline)) {
                TLine* line = new TLine(-0.5, baseline, 6.5, baseline);
                line->SetLineStyle(2);
                line->SetLineWidth(2);
                line->SetLineColor(kGray + 2);
                line->Draw();
                owned.push_back(line);
            }
        }

        TLatex* period_label = new TLatex();
        period_label->SetNDC(true);
        period_label->SetTextAlign(33);
        period_label->SetTextFont(62);
        period_label->SetTextSize(0.040);
        period_label->DrawLatex(0.94, 0.88, period.c_str());
        owned.push_back(period_label);

        if (period == "Sp19 Inb" && hide_sp19_inb_from_all_period_plots) {
            TLatex lat;
            lat.SetNDC(true);
            lat.SetTextSize(0.030);
            lat.SetTextColor(kRed + 2);
            lat.DrawLatex(0.13, 0.88, "Sp19 Inb efficiency copied from Fa18 Inb");
        }
    }

    c.SaveAs(out_path.c_str());
    for (TObject* obj : owned) {
        delete obj;
    }
}

static void draw_photon_region_slope_canvas(const std::string& out_path,
                                            const PhotonRegionResults& results,
                                            bool data_sample,
                                            bool hide_sp19_inb_from_all_period_plots) {
    mkdir_p(out_path.substr(0, out_path.find_last_of('/')));

    TCanvas c(data_sample ? "c_photon_region_data_slope" : "c_photon_region_mc_slope",
              data_sample ? "DATA current slope by photon region" : "MC current slope by photon region",
              1800, 1200);
    c.Divide(3, 2, 0.002, 0.002);

    const std::vector<std::string> panel_order = {
        "Sp18 Inb", "Sp18 Out", "", "Fa18 Inb", "Fa18 Out", "Sp19 Inb"
    };

    std::vector<TObject*> owned;

    for (int ipad = 0; ipad < 6; ++ipad) {
        c.cd(ipad + 1);
        gPad->SetGridy();
        const std::string period = panel_order[ipad];

        if (period.empty()) {
            TLatex lat;
            lat.SetNDC();
            lat.SetTextSize(0.055);
            lat.DrawLatex(0.10, 0.82, data_sample ? "DATA slope" : "MC slope");
            lat.SetTextSize(0.045);
            lat.DrawLatex(0.10, 0.68, "100 m / b  [% / nA]");
            lat.DrawLatex(0.10, 0.54, "Diagnostic FT + FD-sector split");
            continue;
        }

        TH1D* frame = new TH1D(("h_photon_region_slope_" + std::to_string(ipad) + (data_sample ? "_d" : "_m")).c_str(),
                               (period + ";Photon region;Current-dependence slope (%/nA)").c_str(),
                               7, -0.5, 6.5);
        frame->SetStats(0);
        frame->SetMinimum(-1.5);
        frame->SetMaximum(0.5);
        for (int ir = 0; ir < 7; ++ir) {
            frame->GetXaxis()->SetBinLabel(ir + 1, PHOTON_REGION_ORDER[ir].label.c_str());
        }
        frame->Draw("AXIS");
        owned.push_back(frame);

        TLine* zero = new TLine(-0.5, 0.0, 6.5, 0.0);
        zero->SetLineStyle(2);
        zero->SetLineColor(kGray + 2);
        zero->Draw();
        owned.push_back(zero);

        TGraphErrors* graph = new TGraphErrors();
        int ip = 0;
        for (int ir = 0; ir < 7; ++ir) {
            auto it = results.find(PHOTON_REGION_ORDER[ir].key);
            if (it == results.end()) continue;
            const PeriodResult* r = find_period_result(it->second, period);
            if (!r) continue;
            const FitResult& fit = data_sample ? r->data_fit : r->mc_fit;
            const double value = fit_percent_slope(fit);
            const double err = fit_percent_slope_err(fit);
            if (!std::isfinite(value)) continue;
            graph->SetPoint(ip, (double)ir, value);
            graph->SetPointError(ip, 0.0, std::isfinite(err) ? err : 0.0);
            ++ip;
        }
        graph->SetMarkerStyle(20);
        graph->SetMarkerSize(1.2);
        graph->SetLineWidth(2);
        graph->Draw("PE SAME");
        owned.push_back(graph);

        TLatex* period_label = new TLatex();
        period_label->SetNDC(true);
        period_label->SetTextAlign(33);
        period_label->SetTextFont(62);
        period_label->SetTextSize(0.040);
        period_label->DrawLatex(0.94, 0.88, period.c_str());
        owned.push_back(period_label);

        if (period == "Sp19 Inb" && hide_sp19_inb_from_all_period_plots) {
            TLatex lat;
            lat.SetNDC(true);
            lat.SetTextSize(0.030);
            lat.SetTextColor(kRed + 2);
            lat.DrawLatex(0.13, 0.88, "Sp19 Inb efficiency copied from Fa18 Inb");
        }
    }

    c.SaveAs(out_path.c_str());
    for (TObject* obj : owned) {
        delete obj;
    }
}




static double ratio_with_error(double a, double sa,
                               double b, double sb,
                               double& sr) {
    sr = std::numeric_limits<double>::quiet_NaN();
    if (!(std::isfinite(a) && a > 0.0 &&
          std::isfinite(b) && b > 0.0 &&
          std::isfinite(sa) && sa >= 0.0 &&
          std::isfinite(sb) && sb >= 0.0)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    const double r = a / b;
    sr = std::fabs(r) * std::sqrt(
        std::pow(sa / a, 2) + std::pow(sb / b, 2));
    return r;
}

static void write_and_draw_ft_fd_current_double_ratio(
    const std::string& output_dir,
    const std::string& channel_token,
    const PhotonRegionResults& results,
    bool use_fa18_for_sp19) {

    auto it_ft = results.find("FT");
    auto it_fd = results.find("FD");
    if (it_ft == results.end() || it_fd == results.end()) return;

    const std::string odir =
        output_dir + "/" + channel_token + "/sector_dependence_diagnostic";
    mkdir_p(odir);

    const std::string csv_path =
        odir + "/ft_fd_relative_current_efficiency_data_over_mc.csv";
    std::ofstream fout(csv_path);
    fout << "period,ft_data_factor,ft_data_factor_stat,ft_mc_factor,ft_mc_factor_stat,"
         << "fd_data_factor,fd_data_factor_stat,fd_mc_factor,fd_mc_factor_stat,"
         << "ft_data_over_mc,ft_data_over_mc_stat,fd_data_over_mc,fd_data_over_mc_stat,"
         << "ft_over_fd_double_ratio,ft_over_fd_double_ratio_stat\n";

    TCanvas c(("c_ft_fd_current_double_ratio_" + channel_token).c_str(), "", 1250, 800);
    configure_production_pad((TPad*)&c, 0.13, 0.035, 0.16, 0.08);

    TH1D frame(("h_ft_fd_current_double_ratio_" + channel_token).c_str(),
               ";Run period;(DATA/MC)_{FT} / (DATA/MC)_{FD}",
               (int)PERIOD_ORDER.size(), 0.5, PERIOD_ORDER.size() + 0.5);
    frame.SetMinimum(0.50);
    frame.SetMaximum(1.50);
    frame.SetStats(0);
    style_production_frame(&frame, 0.050, 0.050, 0.040);
    for (int i = 0; i < (int)PERIOD_ORDER.size(); ++i)
        frame.GetXaxis()->SetBinLabel(i + 1, PERIOD_ORDER[i].c_str());
    frame.Draw("AXIS");

    TLine unity(0.5, 1.0, PERIOD_ORDER.size() + 0.5, 1.0);
    unity.SetLineStyle(2);
    unity.SetLineWidth(2);
    unity.SetLineColor(kGray + 1);
    unity.Draw("SAME");

    TGraphErrors graph;
    graph.SetMarkerStyle(20);
    graph.SetMarkerSize(1.35);
    graph.SetLineWidth(2);

    int ip = 0;
    for (int i = 0; i < (int)PERIOD_ORDER.size(); ++i) {
        const std::string& period = PERIOD_ORDER[i];
        const PeriodResult* ft = find_period_result(it_ft->second, period);
        const PeriodResult* fd = find_period_result(it_fd->second, period);
        if (!ft || !fd) continue;

        double s_ft_dm = std::numeric_limits<double>::quiet_NaN();
        double s_fd_dm = std::numeric_limits<double>::quiet_NaN();
        const double ft_dm = ratio_with_error(
            ft->data_factor, ft->data_factor_err,
            ft->mc_factor, ft->mc_factor_err, s_ft_dm);
        const double fd_dm = ratio_with_error(
            fd->data_factor, fd->data_factor_err,
            fd->mc_factor, fd->mc_factor_err, s_fd_dm);

        double s_double = std::numeric_limits<double>::quiet_NaN();
        const double double_ratio =
            ratio_with_error(ft_dm, s_ft_dm, fd_dm, s_fd_dm, s_double);

        fout << period << ","
             << ft->data_factor << "," << ft->data_factor_err << ","
             << ft->mc_factor << "," << ft->mc_factor_err << ","
             << fd->data_factor << "," << fd->data_factor_err << ","
             << fd->mc_factor << "," << fd->mc_factor_err << ","
             << ft_dm << "," << s_ft_dm << ","
             << fd_dm << "," << s_fd_dm << ","
             << double_ratio << "," << s_double << "\n";

        if (std::isfinite(double_ratio) && std::isfinite(s_double)) {
            graph.SetPoint(ip, i + 1.0, double_ratio);
            graph.SetPointError(ip, 0.0, s_double);
            ++ip;
        }
    }

    graph.Draw("PE SAME");

    TLatex note;
    note.SetNDC();
    note.SetTextFont(42);
    note.SetTextSize(0.031);
    note.DrawLatex(0.15, 0.89,
                   "unity: FT and combined FD current response are equally modeled");
    if (use_fa18_for_sp19) {
        note.SetTextColor(kRed + 2);
        note.SetTextSize(0.028);
        note.DrawLatex(0.15, 0.84,
                       "Sp19 Inb DATA efficiency copied from Fa18 Inb");
    }

    const std::string png_path =
        odir + "/ft_fd_relative_current_efficiency_data_over_mc.png";
    c.SaveAs(png_path.c_str());

    const std::string note_dir = output_dir + "/analysis_note";
    mkdir_p(note_dir);
    gSystem->CopyFile(
        png_path.c_str(),
        (note_dir + "/" + channel_token +
         "_ft_fd_relative_current_efficiency_data_over_mc.png").c_str(),
        true);
}

static PhotonRegionResults run_photon_region_current_diagnostic(
    const ChannelConfig& cfg,
    const std::map<std::string, TTree*>& data_trees,
    const std::map<std::string, TTree*>& gen_trees,
    const std::map<std::string, TTree*>& rec_trees,
    const std::unordered_map<int, ChargeEntry>& charge_map,
    const TopoCutMap& data_cuts,
    const TopoCutMap& mc_cuts,
    const std::string& output_dir,
    int max_workers,
    bool use_second_column_charge_for_all_unpolarized,
    bool use_columns_3_to_5_charge_sum_scaled_for_fa18_sp19_unpolarized,
    double columns_3_to_5_charge_sum_scale,
    bool use_fa18_inb_current_efficiency_for_sp19_inb,
    bool process_mc,
    const std::vector<PeriodResult>& integrated_results) {

    using Clock = std::chrono::steady_clock;
    const auto diagnostic_t0 = Clock::now();
    auto elapsed_s = [](const Clock::time_point& t0) {
        return std::chrono::duration<double>(Clock::now() - t0).count();
    };

    std::cout << "[current_dependence] Starting photon-region diagnostic for "
              << cfg.csv_channel << ": FT + FD sectors 1--6; "
              << "regional MC scan=" << (process_mc ? "enabled" : "skipped")
              << "." << std::endl;
#ifdef _OPENMP
    std::cout << "[current_dependence] Photon-region diagnostic OpenMP: enabled; requested worker cap="
              << max_workers << ", runtime max threads=" << omp_get_max_threads() << "."
              << std::endl;
#else
    std::cout << "[current_dependence] WARNING: photon-region diagnostic was compiled WITHOUT OpenMP; "
              << "all worker loops will execute serially despite max_workers="
              << max_workers << "." << std::endl;
#endif

    std::map<std::string, std::map<std::string, DataAgg>> data_aggs_by_region_period;
    std::vector<std::pair<std::string, TTree*>> data_items;

    for (const auto& kv : data_trees) {
        PeriodTags tags = parse_period_from_key(kv.first);
        if (tags.display != "Fa18 Inb Supp") {
            data_items.push_back(kv);
        }
    }

    std::mutex data_mutex;
    std::mutex progress_mutex;
    std::atomic<int> data_done{0};
    int nth = std::max(1, std::min(7, max_workers));
    nth = std::min(nth, std::max(1, (int)data_items.size()));

    std::cout << "[current_dependence] Photon-region DATA: " << data_items.size()
              << " tree(s), " << nth << " worker(s)." << std::endl;

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(nth)
#endif
    for (int i = 0; i < (int)data_items.size(); ++i) {
        const auto& item = data_items[i];
        const auto tree_t0 = Clock::now();
        const Long64_t nentries = item.second ? item.second->GetEntries() : 0;

        {
            std::lock_guard<std::mutex> lock(progress_mutex);
            std::cout << "[current_dependence] Photon-region DATA start "
                      << (i + 1) << "/" << data_items.size()
                      << " key=" << item.first
                      << " entries=" << nentries << std::endl;
        }

        const PhotonRegionDataAggMap local = process_data_tree_by_photon_region(
            cfg,
            item.first,
            item.second,
            charge_map,
            data_cuts,
            use_second_column_charge_for_all_unpolarized,
            use_columns_3_to_5_charge_sum_scaled_for_fa18_sp19_unpolarized,
            columns_3_to_5_charge_sum_scale);

        {
            std::lock_guard<std::mutex> lock(data_mutex);
            for (const PhotonRegionSpec& spec : PHOTON_REGION_ORDER) {
                auto it = local.find(spec.key);
                if (it != local.end()) {
                    data_aggs_by_region_period[spec.key][it->second.period] = it->second;
                }
            }
        }

        const int finished = ++data_done;
        {
            std::lock_guard<std::mutex> lock(progress_mutex);
            std::cout << "[current_dependence] Photon-region DATA done "
                      << finished << "/" << data_items.size()
                      << " key=" << item.first
                      << " entries=" << nentries
                      << " elapsed=" << std::fixed << std::setprecision(1)
                      << elapsed_s(tree_t0) << " s" << std::defaultfloat
                      << std::endl;
        }
    }

    std::cout << "[current_dependence] Photon-region DATA phase complete in "
              << std::fixed << std::setprecision(1) << elapsed_s(diagnostic_t0)
              << " s." << std::defaultfloat << std::setprecision(6) << std::endl;

    // MC counts are stored by region and by period/current.  For ep->eppi0
    // this scan is intentionally skipped: the regional pi0 MC factors are
    // constructed downstream from the measured pi0 DATA response and the
    // DVCS regional DATA/MC ratio, so scanning pi0 generated/reconstructed MC
    // here would be expensive work whose result is immediately overwritten.
    std::map<std::string, std::map<std::string, McAgg>> mc_by_region_period_current;

    if (!process_mc) {
        std::cout << "[current_dependence] Photon-region " << cfg.csv_channel
                  << " MC phase skipped by construction; regional MC factors "
                  << "will be supplied by the downstream DVCS DATA/MC transfer."
                  << std::endl;
    } else {
        std::vector<std::pair<std::string, TTree*>> gen_items;
        gen_items.reserve(gen_trees.size());
        for (const auto& kv : gen_trees) {
            const PeriodTags tags = parse_period_from_key(kv.first);
            if (tags.display != "Fa18 Inb Supp") gen_items.push_back(kv);
        }

        std::mutex mc_mutex;
        std::atomic<int> gen_done{0};
        nth = std::max(1, std::min(7, max_workers));
        nth = std::min(nth, std::max(1, (int)gen_items.size()));

        std::cout << "[current_dependence] Photon-region generated MC: "
                  << gen_items.size() << " tree(s), " << nth
                  << " worker(s)." << std::endl;
        const auto gen_t0 = Clock::now();

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(nth)
#endif
        for (int i = 0; i < (int)gen_items.size(); ++i) {
            const auto& item = gen_items[i];
            const PeriodTags tags = parse_period_from_key(item.first);
            const int current = parse_current_from_key(item.first);
            const auto tree_t0 = Clock::now();
            const Long64_t nentries = item.second ? item.second->GetEntries() : 0;

            {
                std::lock_guard<std::mutex> lock(progress_mutex);
                std::cout << "[current_dependence] Photon-region GEN start "
                          << (i + 1) << "/" << gen_items.size()
                          << " key=" << item.first
                          << " entries=" << nentries << std::endl;
            }

            const PhotonRegionMcCountMap counts =
                count_generated_tree_by_photon_region(item.first, item.second);

            {
                std::lock_guard<std::mutex> lock(mc_mutex);
                for (const PhotonRegionSpec& spec : PHOTON_REGION_ORDER) {
                    std::ostringstream k;
                    k << tags.display << "_" << current;
                    McAgg& agg = mc_by_region_period_current[spec.key][k.str()];
                    agg.period = tags.display;
                    agg.current_nA = current;
                    auto it = counts.find(spec.key);
                    if (it != counts.end()) agg.n_gen += it->second;
                }
            }

            const int finished = ++gen_done;
            {
                std::lock_guard<std::mutex> lock(progress_mutex);
                std::cout << "[current_dependence] Photon-region GEN done "
                          << finished << "/" << gen_items.size()
                          << " key=" << item.first
                          << " entries=" << nentries
                          << " elapsed=" << std::fixed << std::setprecision(1)
                          << elapsed_s(tree_t0) << " s" << std::defaultfloat
                          << std::endl;
            }
        }

        std::cout << "[current_dependence] Photon-region generated-MC phase complete in "
                  << std::fixed << std::setprecision(1) << elapsed_s(gen_t0)
                  << " s." << std::defaultfloat << std::setprecision(6) << std::endl;

        std::vector<std::pair<std::string, TTree*>> rec_items;
        rec_items.reserve(rec_trees.size());
        for (const auto& kv : rec_trees) {
            const PeriodTags tags = parse_period_from_key(kv.first);
            if (tags.display != "Fa18 Inb Supp") rec_items.push_back(kv);
        }

        std::atomic<int> rec_done{0};
        nth = std::max(1, std::min(7, max_workers));
        nth = std::min(nth, std::max(1, (int)rec_items.size()));

        std::cout << "[current_dependence] Photon-region reconstructed MC: "
                  << rec_items.size() << " tree(s), " << nth
                  << " worker(s)." << std::endl;
        const auto rec_t0 = Clock::now();

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(nth)
#endif
        for (int i = 0; i < (int)rec_items.size(); ++i) {
            const auto& item = rec_items[i];
            const PeriodTags tags = parse_period_from_key(item.first);
            const int current = parse_current_from_key(item.first);
            const auto tree_t0 = Clock::now();
            const Long64_t nentries = item.second ? item.second->GetEntries() : 0;

            {
                std::lock_guard<std::mutex> lock(progress_mutex);
                std::cout << "[current_dependence] Photon-region REC start "
                          << (i + 1) << "/" << rec_items.size()
                          << " key=" << item.first
                          << " entries=" << nentries << std::endl;
            }

            const PhotonRegionMcCountMap counts =
                count_reconstructed_tree_by_photon_region(
                    cfg, item.first, item.second, mc_cuts);

            {
                std::lock_guard<std::mutex> lock(mc_mutex);
                for (const PhotonRegionSpec& spec : PHOTON_REGION_ORDER) {
                    std::ostringstream k;
                    k << tags.display << "_" << current;
                    McAgg& agg = mc_by_region_period_current[spec.key][k.str()];
                    agg.period = tags.display;
                    agg.current_nA = current;
                    auto it = counts.find(spec.key);
                    if (it != counts.end()) agg.n_rec += it->second;
                }
            }

            const int finished = ++rec_done;
            {
                std::lock_guard<std::mutex> lock(progress_mutex);
                std::cout << "[current_dependence] Photon-region REC done "
                          << finished << "/" << rec_items.size()
                          << " key=" << item.first
                          << " entries=" << nentries
                          << " elapsed=" << std::fixed << std::setprecision(1)
                          << elapsed_s(tree_t0) << " s" << std::defaultfloat
                          << std::endl;
            }
        }

        std::cout << "[current_dependence] Photon-region reconstructed-MC phase complete in "
                  << std::fixed << std::setprecision(1) << elapsed_s(rec_t0)
                  << " s." << std::defaultfloat << std::setprecision(6) << std::endl;
    }

    PhotonRegionResults results;

    for (const PhotonRegionSpec& spec : PHOTON_REGION_ORDER) {
        std::vector<McAgg> mc_aggs;
        for (const auto& kv : mc_by_region_period_current[spec.key]) {
            mc_aggs.push_back(kv.second);
        }

        std::vector<PeriodResult> region_rows;
        for (const std::string& period : PERIOD_ORDER) {
            PeriodResult r;
            r.period = period;

            auto it_data = data_aggs_by_region_period[spec.key].find(period);
            if (it_data != data_aggs_by_region_period[spec.key].end()) {
                r.data_points = data_points_from_agg(it_data->second);
                r.data_fit = fit_points(r.data_points);
                r.data_factor = weighted_data_rel(r.data_points, r.data_fit);
                r.data_factor_err = weighted_data_rel_err(r.data_points, r.data_fit);
            }

            if (process_mc) {
                r.mc_points = mc_points_from_aggs(mc_aggs, period);
                r.mc_fit = fit_points(r.mc_points);
                const int ref = reference_current_nA(period);
                r.mc_factor = rel_at_current((double)ref, r.mc_fit);
                r.mc_factor_err = rel_err_at_current((double)ref, r.mc_fit);
            }

            region_rows.push_back(r);
        }

        if (use_fa18_inb_current_efficiency_for_sp19_inb) {
            replace_sp19_inb_factors_with_fa18_inb(
                region_rows,
                cfg.csv_channel + " photon-region " + spec.label);
        }

        results[spec.key] = std::move(region_rows);
    }

    // Direct combined-FD response.  This is not an average of the six fitted
    // sector factors: the S1--S6 events are merged first and then the DATA and
    // MC current responses are refit.  It is used only for the FT-vs-FD
    // relative-modeling diagnostic below.
    std::vector<PeriodResult> fd_rows;
    for (const std::string& period : PERIOD_ORDER) {
        PeriodResult fd;
        fd.period = period;

        DataAgg fd_data;
        fd_data.period = period;
        for (int is = 1; is <= 6; ++is) {
            const std::string skey = "S" + std::to_string(is);
            auto ir = data_aggs_by_region_period.find(skey);
            if (ir == data_aggs_by_region_period.end()) continue;
            auto iper = ir->second.find(period);
            if (iper != ir->second.end()) merge_data_agg(fd_data, iper->second);
        }
        fd.data_points = data_points_from_agg(fd_data);
        fd.data_fit = fit_points(fd.data_points);
        fd.data_factor = weighted_data_rel(fd.data_points, fd.data_fit);
        fd.data_factor_err = weighted_data_rel_err(fd.data_points, fd.data_fit);

        if (process_mc) {
            std::map<int, McAgg> combined_by_current;
            for (int is = 1; is <= 6; ++is) {
                const std::string skey = "S" + std::to_string(is);
                auto ir = mc_by_region_period_current.find(skey);
                if (ir == mc_by_region_period_current.end()) continue;
                for (const auto& kv : ir->second) {
                    const McAgg& src = kv.second;
                    if (src.period != period) continue;
                    McAgg& dst = combined_by_current[src.current_nA];
                    dst.period = period;
                    dst.current_nA = src.current_nA;
                    dst.n_gen += src.n_gen;
                    dst.n_rec += src.n_rec;
                }
            }
            std::vector<McAgg> combined;
            for (const auto& kv : combined_by_current) combined.push_back(kv.second);
            fd.mc_points = mc_points_from_aggs(combined, period);
            fd.mc_fit = fit_points(fd.mc_points);
            const int ref = reference_current_nA(period);
            fd.mc_factor = rel_at_current((double)ref, fd.mc_fit);
            fd.mc_factor_err = rel_err_at_current((double)ref, fd.mc_fit);
        }

        fd_rows.push_back(fd);
    }

    if (use_fa18_inb_current_efficiency_for_sp19_inb) {
        replace_sp19_inb_factors_with_fa18_inb(
            fd_rows, cfg.csv_channel + " combined FD");
    }
    results["FD"] = fd_rows;

    const std::string odir = output_dir + "/" + cfg.output_token + "/sector_dependence_diagnostic";
    mkdir_p(odir);

    write_photon_region_summary_csv(
        odir + "/photon_region_current_dependence_summary.csv",
        results,
        integrated_results);

    if (process_mc) {
        write_and_draw_ft_fd_current_double_ratio(
            output_dir,
            cfg.output_token,
            results,
            use_fa18_inb_current_efficiency_for_sp19_inb);
    }

    draw_photon_region_response_canvas(
        odir + "/data_current_dependence_by_photon_region.png",
        results,
        true,
        use_fa18_inb_current_efficiency_for_sp19_inb);

    if (process_mc) {
        draw_photon_region_response_canvas(
            odir + "/mc_current_dependence_by_photon_region.png",
            results,
            false,
            use_fa18_inb_current_efficiency_for_sp19_inb);
    }

    draw_photon_region_factor_canvas(
        odir + "/data_current_efficiency_factor_by_photon_region.png",
        results,
        integrated_results,
        true,
        use_fa18_inb_current_efficiency_for_sp19_inb);

    if (process_mc) {
        draw_photon_region_factor_canvas(
            odir + "/mc_current_efficiency_factor_by_photon_region.png",
            results,
            integrated_results,
            false,
            use_fa18_inb_current_efficiency_for_sp19_inb);
    }

    draw_photon_region_slope_canvas(
        odir + "/data_current_slope_by_photon_region.png",
        results,
        true,
        use_fa18_inb_current_efficiency_for_sp19_inb);

    if (process_mc) {
        draw_photon_region_slope_canvas(
            odir + "/mc_current_slope_by_photon_region.png",
            results,
            false,
            use_fa18_inb_current_efficiency_for_sp19_inb);
    }

    std::cout << "[current_dependence] Photon-region diagnostic complete for "
              << cfg.csv_channel << " in "
              << std::fixed << std::setprecision(1) << elapsed_s(diagnostic_t0)
              << " s. Outputs: " << odir << std::defaultfloat << std::setprecision(6) << std::endl;

    return results;
}


static double relative_slope_per_nA(const FitResult& f) {
    if (!std::isfinite(f.m) || !std::isfinite(f.b) || f.b == 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return f.m / f.b;
}

static double relative_slope_per_nA_err(const FitResult& f) {
    if (!std::isfinite(f.m) || !std::isfinite(f.b) ||
        !std::isfinite(f.sm) || !std::isfinite(f.sb) ||
        !std::isfinite(f.cov_mb) || f.b == 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    const double dm = 1.0 / f.b;
    const double db = -f.m / (f.b * f.b);
    double var = dm * dm * f.sm * f.sm + db * db * f.sb * f.sb +
                 2.0 * dm * db * f.cov_mb;
    if (var < 0.0 && std::fabs(var) < 1.0e-15) var = 0.0;
    return (var >= 0.0) ? std::sqrt(var) : std::numeric_limits<double>::quiet_NaN();
}

static void apply_eppi0_mc_factor_from_dvcs_ratio_regions(
    PhotonRegionResults& eppi0_regions,
    const PhotonRegionResults& dvcs_regions) {

    for (const PhotonRegionSpec& spec : PHOTON_REGION_ORDER) {
        auto ie = eppi0_regions.find(spec.key);
        auto id = dvcs_regions.find(spec.key);
        if (ie == eppi0_regions.end() || id == dvcs_regions.end()) continue;

        for (PeriodResult& er : ie->second) {
            const PeriodResult* dr = find_period_result(id->second, er.period);
            if (!dr) continue;
            if (!(std::isfinite(er.data_factor) && er.data_factor > 0.0 &&
                  std::isfinite(dr->data_factor) && dr->data_factor > 0.0 &&
                  std::isfinite(dr->mc_factor) && dr->mc_factor > 0.0)) {
                continue;
            }

            const double dvcs_data_over_mc = dr->data_factor / dr->mc_factor;
            er.mc_factor = er.data_factor / dvcs_data_over_mc;

            const double re = (std::isfinite(er.data_factor_err) && er.data_factor_err >= 0.0)
                                ? er.data_factor_err / er.data_factor : 0.0;
            const double rd = (std::isfinite(dr->data_factor_err) && dr->data_factor_err >= 0.0)
                                ? dr->data_factor_err / dr->data_factor : 0.0;
            const double rm = (std::isfinite(dr->mc_factor_err) && dr->mc_factor_err >= 0.0)
                                ? dr->mc_factor_err / dr->mc_factor : 0.0;
            er.mc_factor_err = er.mc_factor * std::sqrt(re * re + rd * rd + rm * rm);
        }
    }
}

static void collect_run_current_map_from_data_trees(
    const std::map<std::string, TTree*>& trees,
    const std::unordered_map<int, ChargeEntry>& charge_map,
    bool use_second_column_charge_for_all_unpolarized,
    bool use_columns_3_to_5_charge_sum_scaled_for_fa18_sp19_unpolarized,
    double columns_3_to_5_charge_sum_scale,
    nlohmann::json& run_map,
    nlohmann::json& excluded_run_map) {

    for (const auto& kv : trees) {
        if (!kv.second) continue;
        const PeriodTags tags = parse_period_from_key(kv.first);
        if (tags.display == "Fa18 Inb Supp") continue;

        TTree* tree = kv.second;
        int runnum = 0;
        if (!tree->GetBranch("runnum")) continue;

        {
            std::lock_guard<std::mutex> lock(g_bind_mutex);
            tree->SetBranchStatus("*", 0);
            tree->SetBranchStatus("runnum", 1);
            tree->SetBranchAddress("runnum", &runnum);
        }

        const Long64_t N = tree->GetEntries();
        std::unordered_set<int> seen;
        for (Long64_t i = 0; i < N; ++i) {
            tree->GetEntry(i);
            if (!seen.insert(runnum).second) continue;

            // Runs with an explicitly present but unusable normalization
            // charge (most commonly a QADB-rejected/zero-charge run) must not
            // be current-corrected. Record them explicitly so downstream
            // event-level accumulation can skip them without conflating them
            // with a genuinely missing current assignment for a valid run.
            const auto charge_it = charge_map.find(runnum);
            if (charge_it != charge_map.end()) {
                double selected_charge = std::numeric_limits<double>::quiet_NaN();
                const bool charge_ok = select_unpolarized_charge_for_period(
                    tags, runnum, charge_it->second,
                    use_second_column_charge_for_all_unpolarized,
                    use_columns_3_to_5_charge_sum_scaled_for_fa18_sp19_unpolarized,
                    columns_3_to_5_charge_sum_scale, selected_charge);
                if (!charge_ok || !(std::isfinite(selected_charge) && selected_charge > 0.0)) {
                    excluded_run_map[tags.display][std::to_string(runnum)] =
                        "nonpositive_or_unusable_charge";
                    continue;
                }
            }

            int current = 0;
            if (resolve_current(tags.internal, runnum, current)) {
                run_map[tags.display][std::to_string(runnum)] = current;
            }
        }
    }
}


struct FtFdMatchedCell {
    double data_ft = 0.0;
    double data_ft_w2 = 0.0;
    double data_fd = 0.0;
    double data_fd_w2 = 0.0;
    double mc_ft = 0.0;
    double mc_ft_w2 = 0.0;
    double mc_fd = 0.0;
    double mc_fd_w2 = 0.0;
};

struct FtFdMatchedSummary {
    bool valid = false;
    double ratio = std::numeric_limits<double>::quiet_NaN();
    double ratio_stat = std::numeric_limits<double>::quiet_NaN();
    int matched_cells = 0;
};

static FtFdMatchedSummary combine_ft_fd_double_ratios(
    const std::vector<FtFdMatchedCell>& cells) {

    double sw = 0.0;
    double swr = 0.0;
    int n = 0;
    for (const FtFdMatchedCell& c : cells) {
        if (!(c.data_ft > 0.0 && c.data_fd > 0.0 &&
              c.mc_ft > 0.0 && c.mc_fd > 0.0)) continue;
        const double r = (c.data_ft * c.mc_fd) / (c.mc_ft * c.data_fd);
        if (!(std::isfinite(r) && r > 0.0)) continue;
        const double relvar =
            c.data_ft_w2 / (c.data_ft * c.data_ft) +
            c.data_fd_w2 / (c.data_fd * c.data_fd) +
            c.mc_ft_w2 / (c.mc_ft * c.mc_ft) +
            c.mc_fd_w2 / (c.mc_fd * c.mc_fd);
        if (!(std::isfinite(relvar) && relvar > 0.0)) continue;
        const double var = r * r * relvar;
        if (!(std::isfinite(var) && var > 0.0)) continue;
        const double w = 1.0 / var;
        sw += w;
        swr += w * r;
        ++n;
    }

    FtFdMatchedSummary out;
    if (sw > 0.0 && n > 0) {
        out.valid = true;
        out.ratio = swr / sw;
        out.ratio_stat = std::sqrt(1.0 / sw);
        out.matched_cells = n;
    }
    return out;
}

static void run_relative_ft_fd_photon_efficiency_diagnostic(
    const ChannelConfig& cfg,
    const std::map<std::string, TTree*>& data_trees,
    const std::map<std::string, TTree*>& rec_trees,
    const std::unordered_map<int, ChargeEntry>& charge_map,
    const TopoCutMap& data_cuts,
    const TopoCutMap& mc_cuts,
    const PhotonRegionResults& region_results,
    const PeriodPooledAngleFitMap& e_theta_models,
    const std::string& output_dir,
    bool use_second_column_charge_for_all_unpolarized,
    bool use_columns_3_to_5_charge_sum_scaled_for_fa18_sp19_unpolarized,
    double columns_3_to_5_charge_sum_scale,
    bool use_sp18_out_e_theta_response_model) {

    if (cfg.channel != Channel::DVCS) return;

    // Match the electron and proton in a common reconstructed phase space.
    // The proton topology is fixed to CD in both samples, so the comparison is
    // CD_FT versus CD_FD and changes only the photon detector region.
    const std::vector<double> e_edges = {5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0};
    const std::vector<double> p_edges = {10.0, 18.0, 24.0, 30.0, 38.0, 46.0, 54.0, 62.0};
    const int ne = (int)e_edges.size() - 1;
    const int np = (int)p_edges.size() - 1;
    const int ncells = ne * np;
    const std::vector<std::string> periods = {
        "Sp18 Inb", "Sp18 Out", "Fa18 Inb", "Fa18 Out"
    };

    // index 0 = all FD sectors combined; 1..6 = S1..S6.
    using SectorCells = std::array<std::vector<FtFdMatchedCell>, 7>;
    std::map<std::string, SectorCells> by_period;
    for (const std::string& period : periods)
        for (auto& v : by_period[period]) v.resize(ncells);

    auto region_index_from_key = [&](const std::string& key) -> int {
        for (int ir = 0; ir < 7; ++ir)
            if (PHOTON_REGION_ORDER[ir].key == key) return ir;
        return -1;
    };

    auto data_weight = [&](const std::string& period,
                           const std::string& region,
                           int current_nA,
                           double theta_e_deg) -> double {
        auto ir = region_results.find(region);
        if (ir == region_results.end()) return std::numeric_limits<double>::quiet_NaN();
        const PeriodResult* r = find_period_result(ir->second, period);
        if (!r) return std::numeric_limits<double>::quiet_NaN();
        double s = relative_slope_per_nA(r->data_fit);

        if (period == "Sp18 Out" && use_sp18_out_e_theta_response_model) {
            auto im = e_theta_models.find(period);
            const int ridx = region_index_from_key(region);
            if (im != e_theta_models.end() && im->second.valid && ridx >= 0) {
                s = im->second.region_intercepts[ridx] +
                    im->second.theta_gradient *
                    (theta_e_deg - im->second.region_centers[ridx]);
            }
        }
        const double response = 1.0 + s * double(current_nA);
        if (!(std::isfinite(response) && response > 0.0))
            return std::numeric_limits<double>::quiet_NaN();
        return 1.0 / response;
    };

    auto mc_weight = [&](const std::string& period,
                         const std::string& region) -> double {
        auto ir = region_results.find(region);
        if (ir == region_results.end()) return std::numeric_limits<double>::quiet_NaN();
        const PeriodResult* r = find_period_result(ir->second, period);
        if (!r || !(std::isfinite(r->mc_factor) && r->mc_factor > 0.0))
            return std::numeric_limits<double>::quiet_NaN();
        return 1.0 / r->mc_factor;
    };

    // DATA: apply the same current response that production counting will use.
    for (const auto& kv : data_trees) {
        const PeriodTags tags = parse_period_from_key(kv.first);
        if (std::find(periods.begin(), periods.end(), tags.display) == periods.end()) continue;
        TTree* tree = kv.second;
        if (!tree) continue;

        Branches b;
        b.bind(tree);
        if (!(b.has_runnum && b.has_detector1 && b.has_detector2 &&
              b.has_e_theta && b.has_p1_theta)) continue;

        std::unordered_map<int, ResolvedRunInfo> run_cache;
        const Long64_t N = tree->GetEntries();
        for (Long64_t i = 0; i < N; ++i) {
            tree->GetEntry(i);
            if (!passes_cone_cut(b)) continue;
            if (!passes_global_dispatch(b, tags)) continue;
            if (!passes_sigma_dispatch(cfg, tags, data_cuts, b)) continue;

            // Same proton detector topology on both sides of the ratio.
            if (b.detector1 != 2) continue;
            if (!(b.detector2 == 0 || b.detector2 == 1)) continue;

            auto ri = run_cache.find(b.runnum);
            if (ri == run_cache.end()) {
                ri = run_cache.emplace(
                    b.runnum,
                    resolve_run_info_cached(
                        tags, b.runnum, charge_map,
                        use_second_column_charge_for_all_unpolarized,
                        use_columns_3_to_5_charge_sum_scaled_for_fa18_sp19_unpolarized,
                        columns_3_to_5_charge_sum_scale)).first;
            }
            if (!ri->second.valid) continue;

            const double et = b.e_theta * RAD2DEG;
            const double pt = b.p1_theta * RAD2DEG;
            const int ie = find_bin_index(e_edges, et);
            const int ip = find_bin_index(p_edges, pt);
            if (ie < 0 || ip < 0) continue;
            const int icell = ie * np + ip;

            std::string region;
            if (!reconstructed_photon_region(b, region)) continue;
            const double w = data_weight(tags.display, region, ri->second.current, et);
            if (!(std::isfinite(w) && w > 0.0)) continue;

            if (b.detector2 == 0) {
                for (int is = 0; is < 7; ++is) {
                    FtFdMatchedCell& c = by_period[tags.display][is][icell];
                    c.data_ft += w;
                    c.data_ft_w2 += w * w;
                }
            } else {
                const int ridx = region_index_from_key(region);
                if (ridx < 1 || ridx > 6) continue;
                FtFdMatchedCell& all = by_period[tags.display][0][icell];
                all.data_fd += w;
                all.data_fd_w2 += w * w;
                FtFdMatchedCell& sec = by_period[tags.display][ridx][icell];
                sec.data_fd += w;
                sec.data_fd_w2 += w * w;
            }
        }
    }

    // MC: use only the production/reference-current reconstructed sample and
    // remove its regional current-overlay response.  The arbitrary MC
    // normalization cancels in the FT/FD double ratio.
    for (const auto& kv : rec_trees) {
        const PeriodTags tags = parse_period_from_key(kv.first);
        if (std::find(periods.begin(), periods.end(), tags.display) == periods.end()) continue;
        if (parse_current_from_key(kv.first) != reference_current_nA(tags.display)) continue;
        TTree* tree = kv.second;
        if (!tree) continue;

        Branches b;
        b.bind(tree);
        if (!(b.has_detector1 && b.has_detector2 && b.has_e_theta && b.has_p1_theta)) continue;
        const Long64_t N = tree->GetEntries();
        for (Long64_t i = 0; i < N; ++i) {
            tree->GetEntry(i);
            if (!passes_cone_cut(b)) continue;
            if (!passes_global_dispatch(b, tags)) continue;
            if (!passes_sigma_dispatch(cfg, tags, mc_cuts, b)) continue;
            if (b.detector1 != 2) continue;
            if (!(b.detector2 == 0 || b.detector2 == 1)) continue;

            const double et = b.e_theta * RAD2DEG;
            const double pt = b.p1_theta * RAD2DEG;
            const int ie = find_bin_index(e_edges, et);
            const int ip = find_bin_index(p_edges, pt);
            if (ie < 0 || ip < 0) continue;
            const int icell = ie * np + ip;

            std::string region;
            if (!reconstructed_photon_region(b, region)) continue;
            const double w = mc_weight(tags.display, region);
            if (!(std::isfinite(w) && w > 0.0)) continue;

            if (b.detector2 == 0) {
                for (int is = 0; is < 7; ++is) {
                    FtFdMatchedCell& c = by_period[tags.display][is][icell];
                    c.mc_ft += w;
                    c.mc_ft_w2 += w * w;
                }
            } else {
                const int ridx = region_index_from_key(region);
                if (ridx < 1 || ridx > 6) continue;
                FtFdMatchedCell& all = by_period[tags.display][0][icell];
                all.mc_fd += w;
                all.mc_fd_w2 += w * w;
                FtFdMatchedCell& sec = by_period[tags.display][ridx][icell];
                sec.mc_fd += w;
                sec.mc_fd_w2 += w * w;
            }
        }
    }

    const std::string odir = output_dir + "/epg/relative_photon_efficiency";
    const std::string note_dir = output_dir + "/analysis_note";
    mkdir_p(odir);
    mkdir_p(note_dir);

    std::ofstream cell_csv(odir + "/ft_fd_matched_cells.csv");
    cell_csv
        << "period,fd_region,e_theta_low,e_theta_high,p_theta_low,p_theta_high,"
        << "data_ft,data_fd,mc_ft,mc_fd,double_ratio,double_ratio_stat\n";

    std::ofstream summary_csv(odir + "/ft_fd_relative_efficiency_summary.csv");
    summary_csv
        << "period,fd_region,relative_ft_over_fd_data_mc,stat,matched_cells\n";

    std::map<std::string, std::array<FtFdMatchedSummary, 7>> summaries;
    for (const std::string& period : periods) {
        for (int is = 0; is < 7; ++is) {
            summaries[period][is] = combine_ft_fd_double_ratios(by_period[period][is]);
            const std::string label = is == 0 ? "FD_all" : PHOTON_REGION_ORDER[is].key;
            const FtFdMatchedSummary& s = summaries[period][is];
            summary_csv << period << "," << label << ","
                        << s.ratio << "," << s.ratio_stat << ","
                        << s.matched_cells << "\n";

            for (int ie = 0; ie < ne; ++ie) {
                for (int ip = 0; ip < np; ++ip) {
                    const FtFdMatchedCell& c = by_period[period][is][ie*np + ip];
                    double r = std::numeric_limits<double>::quiet_NaN();
                    double se = std::numeric_limits<double>::quiet_NaN();
                    if (c.data_ft > 0.0 && c.data_fd > 0.0 &&
                        c.mc_ft > 0.0 && c.mc_fd > 0.0) {
                        r = (c.data_ft * c.mc_fd) / (c.mc_ft * c.data_fd);
                        const double relvar =
                            c.data_ft_w2/(c.data_ft*c.data_ft) +
                            c.data_fd_w2/(c.data_fd*c.data_fd) +
                            c.mc_ft_w2/(c.mc_ft*c.mc_ft) +
                            c.mc_fd_w2/(c.mc_fd*c.mc_fd);
                        if (std::isfinite(relvar) && relvar >= 0.0)
                            se = std::fabs(r) * std::sqrt(relvar);
                    }
                    cell_csv << period << "," << label << ","
                             << e_edges[ie] << "," << e_edges[ie+1] << ","
                             << p_edges[ip] << "," << p_edges[ip+1] << ","
                             << c.data_ft << "," << c.data_fd << ","
                             << c.mc_ft << "," << c.mc_fd << ","
                             << r << "," << se << "\n";
                }
            }
        }
    }

    // Overall period summary: CD_FT relative to all CD_FD sectors.
    TCanvas coverall("c_ft_fd_relative_photon_efficiency", "", 1100, 760);
    configure_production_pad((TPad*)&coverall, 0.14, 0.04, 0.16, 0.08);
    auto oframe = std::make_unique<TH1D>(
        "h_ft_fd_relative_photon_efficiency",
        ";Run period;(DATA/MC)_{CD,FT}/(DATA/MC)_{CD,FD}",
        (int)periods.size(), 0.5, periods.size()+0.5);
    oframe->SetDirectory(nullptr);
    for (size_t ip = 0; ip < periods.size(); ++ip)
        oframe->GetXaxis()->SetBinLabel((int)ip+1, periods[ip].c_str());
    oframe->SetMinimum(0.50);
    oframe->SetMaximum(1.50);
    style_production_frame(oframe.get(), 0.052, 0.052, 0.041);
    oframe->Draw();
    auto unity = std::make_unique<TLine>(0.5, 1.0, periods.size()+0.5, 1.0);
    unity->SetLineStyle(2);
    unity->SetLineWidth(2);
    unity->SetLineColor(kGray+1);
    unity->Draw("SAME");
    auto og = std::make_unique<TGraphErrors>();
    int nog = 0;
    for (size_t ip = 0; ip < periods.size(); ++ip) {
        const auto& s = summaries[periods[ip]][0];
        if (!s.valid) continue;
        og->SetPoint(nog, ip+1.0, s.ratio);
        og->SetPointError(nog, 0.0, s.ratio_stat);
        ++nog;
    }
    og->SetMarkerStyle(20);
    og->SetMarkerSize(1.3);
    og->SetLineWidth(2);
    og->Draw("P SAME");
    TLatex overall_note;
    overall_note.SetNDC();
    overall_note.SetTextFont(42);
    overall_note.SetTextSize(0.032);
    overall_note.DrawLatex(0.16, 0.88, "Matched in (#theta_{e}, #theta_{p}); proton fixed to CD");
    overall_note.DrawLatex(0.16, 0.83, "unity = equal FT/FD photon-efficiency modeling");
    coverall.Modified();
    coverall.Update();
    const std::string overall_path = odir + "/ft_fd_relative_efficiency_summary.png";
    coverall.SaveAs(overall_path.c_str());
    gSystem->CopyFile(
        overall_path.c_str(),
        (note_dir + "/epg_relative_ft_fd_photon_efficiency.png").c_str(), true);

    // Sector-resolved FD comparison.
    TCanvas csector("c_ft_fd_relative_photon_efficiency_sector", "", 1800, 1200);
    csector.Divide(2,2,0.002,0.002);
    std::vector<std::unique_ptr<TObject>> owned;
    for (size_t ip = 0; ip < periods.size(); ++ip) {
        csector.cd((int)ip+1);
        configure_production_pad((TPad*)gPad, 0.145, 0.035, 0.16, 0.075);
        auto frame = std::make_unique<TH1D>(
            ("h_ftfd_sector_"+std::to_string(ip)).c_str(),
            ";FD photon region;(DATA/MC)_{CD,FT}/(DATA/MC)_{CD,FD}",
            7, 0.5, 7.5);
        frame->SetDirectory(nullptr);
        frame->GetXaxis()->SetBinLabel(1, "FD all");
        for (int is=1; is<7; ++is)
            frame->GetXaxis()->SetBinLabel(is+1, PHOTON_REGION_ORDER[is].label.c_str());
        frame->SetMinimum(0.50);
        frame->SetMaximum(1.50);
        style_production_frame(frame.get(), 0.050, 0.050, 0.040);
        frame->Draw();
        owned.emplace_back(std::move(frame));
        auto line = std::make_unique<TLine>(0.5,1.0,7.5,1.0);
        line->SetLineStyle(2); line->SetLineWidth(2); line->SetLineColor(kGray+1);
        line->Draw("SAME");
        owned.emplace_back(std::move(line));
        auto g = std::make_unique<TGraphErrors>();
        int n=0;
        for (int is=0; is<7; ++is) {
            const auto& s=summaries[periods[ip]][is];
            if (!s.valid) continue;
            g->SetPoint(n,is+1.0,s.ratio);
            g->SetPointError(n,0.0,s.ratio_stat);
            ++n;
        }
        g->SetMarkerStyle(20); g->SetMarkerSize(1.25); g->SetLineWidth(2);
        g->Draw("P SAME");
        owned.emplace_back(std::move(g));
        owned.emplace_back(make_period_label(periods[ip],0.94,0.86,33,0.052));
    }
    csector.cd(0);
    auto title=std::make_unique<TLatex>();
    title->SetNDC(); title->SetTextAlign(22); title->SetTextFont(62); title->SetTextSize(0.030);
    title->DrawLatex(0.5,0.970,"Relative FT/FD photon-efficiency modeling after matched charged-particle kinematics");
    owned.emplace_back(std::move(title));
    csector.Modified(); csector.Update();
    const std::string sector_path=odir+"/ft_fd_relative_efficiency_by_fd_sector.png";
    csector.SaveAs(sector_path.c_str());
    gSystem->CopyFile(
        sector_path.c_str(),
        (note_dir+"/epg_relative_ft_fd_photon_efficiency_by_sector.png").c_str(),true);

    std::cout << "[current_dependence] Relative FT/FD photon-efficiency diagnostic written under "
              << odir << std::endl;
}

static void write_current_response_model_json(
    const std::string& path,
    const PhotonRegionResults& dvcs_regions,
    const PhotonRegionResults& eppi0_regions,
    const std::vector<PeriodResult>& dvcs_integrated,
    const std::vector<PeriodResult>& eppi0_integrated,
    const std::map<std::string, TTree*>& dvcs_data_trees,
    const std::map<std::string, TTree*>& eppi0_data_trees,
    const std::unordered_map<int, ChargeEntry>& charge_map,
    bool use_second_column_charge_for_all_unpolarized,
    bool use_columns_3_to_5_charge_sum_scaled_for_fa18_sp19_unpolarized,
    double columns_3_to_5_charge_sum_scale,
    bool use_nobkg_dvcs_mc_counts,
    const PeriodPooledAngleFitMap& dvcs_region_e_theta_models,
    bool use_sp18_out_e_theta_response_model) {

    nlohmann::json j;
    j["schema_version"] = 2;
    j["description"] = "Regional beam-current reconstruction response used for event-level yield correction, with an optional centered polar-angle term";
    j["data_response_form"] = "R(I,theta)=1+I*[s_region + a*(theta-theta_center)] when angular_model is present; otherwise R(I)=1+s_region*I; event weight=1/R";
    j["mc_response_form"] = "reference_factor at production-current MC; event weight=1/reference_factor";
    j["dvcs_nobkg_mc_forced_unity"] = use_nobkg_dvcs_mc_counts;

    auto write_channel = [&](const std::string& channel,
                             const PhotonRegionResults& results,
                             const std::vector<PeriodResult>& integrated) {
        for (const PhotonRegionSpec& spec : PHOTON_REGION_ORDER) {
            auto ir = results.find(spec.key);
            if (ir == results.end()) continue;
            for (const PeriodResult& r : ir->second) {
                const PeriodResult* ri = find_period_result(integrated, r.period);

                double srel = relative_slope_per_nA(r.data_fit);
                double srel_err = relative_slope_per_nA_err(r.data_fit);
                std::string data_source = "regional";
                if (!(std::isfinite(srel) && std::isfinite(srel_err) && srel_err >= 0.0) && ri) {
                    srel = relative_slope_per_nA(ri->data_fit);
                    srel_err = relative_slope_per_nA_err(ri->data_fit);
                    data_source = "integrated_fallback";
                }
                if (!(std::isfinite(srel) && std::isfinite(srel_err) && srel_err >= 0.0)) {
                    fatal("[current_dependence] FATAL: no valid DATA current response for " + channel + " " + r.period + " " + spec.key);
                }

                // Production angular model: only Sp18 Out ep->epg uses the
                // robust common theta_e gradient.  Other periods/channels stay
                // regional-only.  This is an explicit analysis choice rather
                // than an automatic best-BIC switch.
                bool use_angular = false;
                PooledRegionThetaSlopeFit angular_fit;
                int angular_region = -1;
                if (channel == "ep->epg" && r.period == "Sp18 Out" &&
                    use_sp18_out_e_theta_response_model) {
                    auto ia = dvcs_region_e_theta_models.find(r.period);
                    for (int ir = 0; ir < 7; ++ir)
                        if (PHOTON_REGION_ORDER[ir].key == spec.key) angular_region = ir;
                    if (ia != dvcs_region_e_theta_models.end() && ia->second.valid &&
                        angular_region >= 0) {
                        angular_fit = ia->second;
                        srel = angular_fit.region_intercepts[angular_region];
                        srel_err = angular_fit.region_intercept_err[angular_region];
                        data_source = "regional_plus_common_linear_e_theta";
                        use_angular = std::isfinite(srel) && std::isfinite(srel_err) &&
                                      std::isfinite(angular_fit.region_centers[angular_region]) &&
                                      std::isfinite(angular_fit.theta_gradient) &&
                                      std::isfinite(angular_fit.theta_gradient_err);
                    }
                    if (!use_angular) {
                        fatal("[current_dependence] FATAL: Sp18 Out theta_e response model was requested but no valid pooled fit is available for region " + spec.key);
                    }
                }

                auto& d = j["data"][channel][r.period][spec.key];
                d["relative_slope_per_nA"] = srel;
                d["relative_slope_stat"] = srel_err;
                d["source"] = data_source;
                if (use_angular) {
                    d["angular_model"]["variable"] = "e_theta";
                    d["angular_model"]["units"] = "deg";
                    d["angular_model"]["center"] = angular_fit.region_centers[angular_region];
                    d["angular_model"]["gradient_per_nA_per_deg"] = angular_fit.theta_gradient;
                    d["angular_model"]["gradient_stat"] = angular_fit.theta_gradient_err;
                    d["angular_model"]["fit_chi2"] = angular_fit.chi2;
                    d["angular_model"]["fit_ndf"] = angular_fit.ndf;
                }
                if (std::isfinite(r.data_factor)) d["integrated_factor"] = r.data_factor;
                if (std::isfinite(r.data_factor_err)) d["integrated_factor_stat"] = r.data_factor_err;

                double mf = r.mc_factor;
                double me = r.mc_factor_err;
                std::string mc_source = "regional";
                if (channel == "ep->epg" && use_nobkg_dvcs_mc_counts) {
                    mf = 1.0; me = 0.0; mc_source = "nobkg_unity";
                } else if (!(std::isfinite(mf) && mf > 0.0 && std::isfinite(me) && me >= 0.0) && ri) {
                    mf = ri->mc_factor; me = ri->mc_factor_err; mc_source = "integrated_fallback";
                }
                if (!(std::isfinite(mf) && mf > 0.0 && std::isfinite(me) && me >= 0.0)) {
                    fatal("[current_dependence] FATAL: no valid MC current response for " + channel + " " + r.period + " " + spec.key);
                }
                auto& m = j["mc"][channel][r.period][spec.key];
                m["reference_current_nA"] = reference_current_nA(r.period);
                m["reference_factor"] = mf;
                m["reference_factor_stat"] = me;
                m["source"] = mc_source;
            }
        }
    };

    write_channel("ep->epg", dvcs_regions, dvcs_integrated);
    write_channel("ep->eppi0", eppi0_regions, eppi0_integrated);

    collect_run_current_map_from_data_trees(
        dvcs_data_trees, charge_map,
        use_second_column_charge_for_all_unpolarized,
        use_columns_3_to_5_charge_sum_scaled_for_fa18_sp19_unpolarized,
        columns_3_to_5_charge_sum_scale,
        j["run_current_nA"], j["excluded_data_runs"]);
    collect_run_current_map_from_data_trees(
        eppi0_data_trees, charge_map,
        use_second_column_charge_for_all_unpolarized,
        use_columns_3_to_5_charge_sum_scaled_for_fa18_sp19_unpolarized,
        columns_3_to_5_charge_sum_scale,
        j["run_current_nA"], j["excluded_data_runs"]);

    const size_t slash_pos = path.find_last_of('/');
    if (slash_pos != std::string::npos) mkdir_p(path.substr(0, slash_pos));
    std::ofstream fout(path);
    if (!fout.is_open()) fatal("[current_dependence] FATAL: cannot write response model JSON: " + path);
    fout << std::setw(2) << j << "\n";
    std::cout << "[current_dependence] Wrote event-level current-response model: " << path << std::endl;
}

static void draw_analysis_note_region_ratio_canvas(
    const std::string& out_path,
    const PhotonRegionResults& results,
    const std::vector<PeriodResult>& integrated_results,
    const std::string& channel_title,
    bool /*hide_sp19*/) {

    const std::vector<std::string> periods = {"Sp18 Inb", "Sp18 Out", "Fa18 Inb", "Fa18 Out"};
    TCanvas c(("c_note_region_ratio_" + channel_title).c_str(), "", 1800, 1200);
    c.Divide(2, 2, 0.002, 0.002);
    std::vector<std::unique_ptr<TObject>> owned;

    for (size_t ip = 0; ip < periods.size(); ++ip) {
        c.cd((int)ip + 1);
        configure_production_pad((TPad*)gPad, 0.145, 0.035, 0.145, 0.075);

        auto frame = std::make_unique<TH1D>(
            ("h_note_region_ratio_" + std::to_string(ip)).c_str(),
            ";Photon region;(f^{DATA}_{r}/f^{MC}_{r})/(f^{DATA}_{int}/f^{MC}_{int})",
            7, 0.5, 7.5);
        frame->SetDirectory(nullptr);
        for (int ir = 0; ir < 7; ++ir)
            frame->GetXaxis()->SetBinLabel(ir + 1, PHOTON_REGION_ORDER[ir].label.c_str());
        frame->SetMinimum(0.65);
        frame->SetMaximum(1.35);
        style_production_frame(frame.get(), 0.050, 0.050, 0.041);
        frame->Draw();
        owned.emplace_back(std::move(frame));

        const PeriodResult* integ = find_period_result(integrated_results, periods[ip]);
        double cint = std::numeric_limits<double>::quiet_NaN();
        if (integ && std::isfinite(integ->data_factor) && integ->data_factor > 0.0 &&
            std::isfinite(integ->mc_factor) && integ->mc_factor > 0.0) {
            cint = integ->data_factor / integ->mc_factor;
        }

        auto g = std::make_unique<TGraphErrors>();
        int n = 0;
        for (int ir = 0; ir < 7; ++ir) {
            auto it = results.find(PHOTON_REGION_ORDER[ir].key);
            if (it == results.end()) continue;
            const PeriodResult* r = find_period_result(it->second, periods[ip]);
            if (!r || !std::isfinite(cint) || cint <= 0.0 ||
                !(std::isfinite(r->data_factor) && r->data_factor > 0.0 &&
                  std::isfinite(r->mc_factor) && r->mc_factor > 0.0)) continue;

            const double cr = r->data_factor / r->mc_factor;
            const double y = cr / cint;
            const double rel2 = std::pow(r->mc_factor_err / r->mc_factor, 2) +
                                std::pow(r->data_factor_err / r->data_factor, 2);
            g->SetPoint(n, ir + 1, y);
            g->SetPointError(n, 0.0, std::fabs(y) * std::sqrt(std::max(0.0, rel2)));
            ++n;
        }
        g->SetMarkerStyle(20);
        g->SetMarkerSize(1.15);
        g->SetLineWidth(2);
        g->Draw("P SAME");
        owned.emplace_back(std::move(g));

        owned.emplace_back(make_period_label(periods[ip], 0.94, 0.86, 33, 0.052));

        auto unity = std::make_unique<TLine>(0.5, 1.0, 7.5, 1.0);
        unity->SetLineStyle(2);
        unity->SetLineWidth(2);
        unity->Draw("SAME");
        owned.emplace_back(std::move(unity));
    }

    c.cd(0);
    auto title = std::make_unique<TLatex>();
    title->SetNDC();
    title->SetTextAlign(22);
    title->SetTextSize(0.030);
    title->DrawLatex(
        0.5, 0.970,
        (channel_title + ": regional DATA/MC current-efficiency ratio relative to integrated treatment").c_str());
    owned.emplace_back(std::move(title));

    c.Modified();
    c.Update();
    c.SaveAs(out_path.c_str());
}

static void draw_analysis_note_region_absolute_canvas(
    const std::string& out_path,
    const PhotonRegionResults& results,
    const std::string& channel_title) {

    const std::vector<std::string> periods = {"Sp18 Inb", "Sp18 Out", "Fa18 Inb", "Fa18 Out"};
    TCanvas c(("c_note_region_absolute_" + channel_title).c_str(), "", 1800, 1200);
    c.Divide(2, 2, 0.002, 0.002);
    std::vector<std::unique_ptr<TObject>> owned;

    // Fixed production scale: an isolated low-statistics point is allowed
    // to clip rather than determining the scale of all four panels.
    const auto absolute_range =
        production_efficiency_plot_range("data_over_mc");
    const double ymin = absolute_range.first;
    const double ymax = absolute_range.second;

    for (size_t ip = 0; ip < periods.size(); ++ip) {
        c.cd((int)ip + 1);
        configure_production_pad((TPad*)gPad, 0.145, 0.035, 0.145, 0.075);

        auto frame = std::make_unique<TH1D>(
            ("h_note_region_absolute_" + std::to_string(ip)).c_str(),
            ";Photon region;f^{DATA}_{r}/f^{MC}_{r}",
            7, 0.5, 7.5);
        frame->SetDirectory(nullptr);
        for (int ir = 0; ir < 7; ++ir)
            frame->GetXaxis()->SetBinLabel(ir + 1, PHOTON_REGION_ORDER[ir].label.c_str());
        frame->SetMinimum(ymin);
        frame->SetMaximum(ymax);
        style_production_frame(frame.get(), 0.050, 0.050, 0.041);
        frame->Draw();
        owned.emplace_back(std::move(frame));

        auto g = std::make_unique<TGraphErrors>();
        int n = 0;
        for (int ir = 0; ir < 7; ++ir) {
            auto it = results.find(PHOTON_REGION_ORDER[ir].key);
            if (it == results.end()) continue;
            const PeriodResult* r = find_period_result(it->second, periods[ip]);
            if (!r || !(std::isfinite(r->data_factor) && r->data_factor > 0.0 &&
                        std::isfinite(r->mc_factor) && r->mc_factor > 0.0)) continue;
            const double cr = r->data_factor / r->mc_factor;
            const double rel2 = std::pow(r->mc_factor_err / r->mc_factor, 2) +
                                std::pow(r->data_factor_err / r->data_factor, 2);
            g->SetPoint(n, ir + 1, cr);
            g->SetPointError(n, 0.0, std::fabs(cr) * std::sqrt(std::max(0.0, rel2)));
            ++n;
        }
        g->SetMarkerStyle(20);
        g->SetMarkerSize(1.15);
        g->SetLineWidth(2);
        g->Draw("P SAME");
        owned.emplace_back(std::move(g));

        owned.emplace_back(make_period_label(periods[ip], 0.94, 0.86, 33, 0.052));
    }

    c.cd(0);
    auto title = std::make_unique<TLatex>();
    title->SetNDC();
    title->SetTextAlign(22);
    title->SetTextSize(0.030);
    title->DrawLatex(
        0.5, 0.970,
        (channel_title + ": nominal regional DATA/MC current-efficiency ratio").c_str());
    owned.emplace_back(std::move(title));

    c.Modified();
    c.Update();
    c.SaveAs(out_path.c_str());
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

        const double dvcs_data_over_mc = dvcs.data_factor / dvcs.mc_factor;
        const double corrected = eppi0.data_factor / dvcs_data_over_mc;

        const double rel_var =
            std::pow(eppi0.data_factor_err / eppi0.data_factor, 2.0) +
            std::pow(dvcs.mc_factor_err / dvcs.mc_factor, 2.0) +
            std::pow(dvcs.data_factor_err / dvcs.data_factor, 2.0);

        const double corrected_err = std::fabs(corrected) * std::sqrt(rel_var);

        eppi0.mc_factor = corrected;
        eppi0.mc_factor_err = corrected_err;

        std::cout << "[current_dependence] Built eppi0 MC current factor for "
                  << eppi0.period
                  << " using eppi0_data / (dvcs_data / dvcs_mc): "
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

#ifdef _OPENMP
        std::cout << "[current_dependence] OpenMP enabled: max_workers="
                  << options.max_workers
                  << ", omp_get_max_threads()=" << omp_get_max_threads()
                  << "." << std::endl;
#else
        std::cout << "[current_dependence] WARNING: executable was built without OpenMP; "
                  << "current-dependence tree loops will be serial." << std::endl;
#endif

        CSV csv;
        load_csv_strict(csv_path, csv);

        const ChannelConfig dvcs = dvcs_config();
        const ChannelConfig eppi0 = eppi0_config();

        PeriodLinearFitMap dvcs_e_theta_data_fits;
        PeriodLinearFitMap eppi0_e_theta_data_fits;
        PeriodPooledAngleFitMap dvcs_region_e_theta_models;

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

        std::cout << "[current_dependence] Loading Python-optimized nominal "
                  << "production exclusivity cuts from "
                  << options.combined_cuts_json << std::endl;

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

        PhotonRegionResults dvcs_region_results;
        PhotonRegionResults eppi0_region_results;

        if (options.enable_photon_region_current_diagnostic) {
            dvcs_region_results = run_photon_region_current_diagnostic(
                dvcs,
                dvcsDataTrees,
                dvcsGenMcTrees,
                dvcsRecMcTrees,
                charge_map,
                data_cuts,
                mc_cuts,
                options.output_dir,
                options.max_workers,
                options.use_second_column_charge_for_all_unpolarized,
                options.use_columns_3_to_5_charge_sum_scaled_for_fa18_sp19_unpolarized,
                options.columns_3_to_5_charge_sum_scale,
                options.use_fa18_inb_current_efficiency_for_sp19_inb,
                true,
                dvcs_results);
        }

        if (options.enable_eppi0_photon_region_current_diagnostic) {
            eppi0_region_results = run_photon_region_current_diagnostic(
                eppi0,
                eppi0DataTrees,
                eppi0GenMcTrees,
                eppi0RecMcTrees,
                charge_map,
                data_cuts,
                mc_cuts,
                options.output_dir,
                options.max_workers,
                options.use_second_column_charge_for_all_unpolarized,
                options.use_columns_3_to_5_charge_sum_scaled_for_fa18_sp19_unpolarized,
                options.columns_3_to_5_charge_sum_scale,
                options.use_fa18_inb_current_efficiency_for_sp19_inb,
                false,
                eppi0_results);
        }

        if (options.enable_region_theta_current_diagnostic) {
            dvcs_region_e_theta_models = run_region_theta_data_diagnostic(
                dvcs, dvcsDataTrees, dvcsGenMcTrees, dvcsRecMcTrees,
                charge_map, data_cuts, mc_cuts, options.output_dir, options.max_workers,
                options.use_second_column_charge_for_all_unpolarized,
                options.use_columns_3_to_5_charge_sum_scaled_for_fa18_sp19_unpolarized,
                options.columns_3_to_5_charge_sum_scale,
                options.use_fa18_inb_current_efficiency_for_sp19_inb,
                true);
            (void)run_region_theta_data_diagnostic(
                eppi0, eppi0DataTrees, eppi0GenMcTrees, eppi0RecMcTrees,
                charge_map, data_cuts, mc_cuts, options.output_dir, options.max_workers,
                options.use_second_column_charge_for_all_unpolarized,
                options.use_columns_3_to_5_charge_sum_scaled_for_fa18_sp19_unpolarized,
                options.columns_3_to_5_charge_sum_scale,
                options.use_fa18_inb_current_efficiency_for_sp19_inb,
                false);
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

        if (!dvcs_region_results.empty() && !eppi0_region_results.empty()) {
            apply_eppi0_mc_factor_from_dvcs_ratio_regions(eppi0_region_results, dvcs_region_results);

            // The direct-pi0 MC scan has only the production-current sample, so
            // overwrite its preliminary unity regional MC factors with the same
            // DVCS DATA/MC transfer prescription used for the integrated pi0
            // correction, now evaluated region by region.
            const std::string eppi0_region_dir = options.output_dir + "/" + eppi0.output_token + "/sector_dependence_diagnostic";
            write_photon_region_summary_csv(eppi0_region_dir + "/photon_region_current_dependence_summary.csv",
                                            eppi0_region_results, eppi0_results);
            draw_photon_region_factor_canvas(eppi0_region_dir + "/mc_current_efficiency_factor_by_photon_region.png",
                                             eppi0_region_results, eppi0_results, false,
                                             options.use_fa18_inb_current_efficiency_for_sp19_inb);
        }

        if (options.use_fa18_inb_current_efficiency_for_sp19_inb) {
            replace_sp19_inb_factors_with_fa18_inb(eppi0_results, eppi0.csv_channel);
        }

        if (options.enable_relative_ft_fd_photon_efficiency_diagnostic &&
            !dvcs_region_results.empty()) {
            run_relative_ft_fd_photon_efficiency_diagnostic(
                dvcs, dvcsDataTrees, dvcsRecMcTrees, charge_map,
                data_cuts, mc_cuts, dvcs_region_results,
                dvcs_region_e_theta_models, options.output_dir,
                options.use_second_column_charge_for_all_unpolarized,
                options.use_columns_3_to_5_charge_sum_scaled_for_fa18_sp19_unpolarized,
                options.columns_3_to_5_charge_sum_scale,
                options.use_sp18_out_e_theta_response_model);
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

        if (!dvcs_region_results.empty() && !eppi0_region_results.empty()) {
            write_current_response_model_json(options.response_model_json,
                                              dvcs_region_results,
                                              eppi0_region_results,
                                              dvcs_results,
                                              eppi0_results,
                                              dvcsDataTrees,
                                              eppi0DataTrees,
                                              charge_map,
                                              options.use_second_column_charge_for_all_unpolarized,
                                              options.use_columns_3_to_5_charge_sum_scaled_for_fa18_sp19_unpolarized,
                                              options.columns_3_to_5_charge_sum_scale,
                                              options.use_nobkg_dvcs_mc_counts,
                                              dvcs_region_e_theta_models,
                                              options.use_sp18_out_e_theta_response_model);

            const std::string note_dir = options.output_dir + "/analysis_note";
            mkdir_p(note_dir);
            draw_all_period_data_mc_canvas(note_dir + "/epg_current_response_data_mc.png",
                                           "ep #rightarrow ep#gamma current response",
                                           dvcs_results,
                                           options.use_fa18_inb_current_efficiency_for_sp19_inb);
            draw_all_period_data_mc_canvas(note_dir + "/eppi0_current_response_data_mc.png",
                                           "ep #rightarrow ep#pi^{0} current response",
                                           eppi0_results,
                                           options.use_fa18_inb_current_efficiency_for_sp19_inb);
            draw_analysis_note_region_ratio_canvas(note_dir + "/epg_regional_current_correction.png",
                                                   dvcs_region_results, dvcs_results,
                                                   "ep #rightarrow ep#gamma",
                                                   options.use_fa18_inb_current_efficiency_for_sp19_inb);
            draw_analysis_note_region_ratio_canvas(note_dir + "/eppi0_regional_current_correction.png",
                                                   eppi0_region_results, eppi0_results,
                                                   "ep #rightarrow ep#pi^{0}",
                                                   options.use_fa18_inb_current_efficiency_for_sp19_inb);
            draw_analysis_note_region_absolute_canvas(
                note_dir + "/epg_regional_current_correction_absolute.png",
                dvcs_region_results, "ep #rightarrow ep#gamma");
            draw_analysis_note_region_absolute_canvas(
                note_dir + "/eppi0_regional_current_correction_absolute.png",
                eppi0_region_results, "ep #rightarrow ep#pi^{0}");
        }

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
        if (options.apply_legacy_binned_current_corrections) {
            apply_all_data_current_corrections(csv,
                                               dvcs,
                                               eppi0,
                                               dvcs_e_theta_data_fits,
                                               eppi0_e_theta_data_fits,
                                               options.use_e_theta_linear_data_current_efficiency);

            // Legacy compatibility mode only. Nominal production fills these
            // corrected columns event-by-event in total_counts.cpp.
            apply_all_mc_current_corrections(csv,
                                         dvcs,
                                         eppi0,
                                         options.use_epg_mc_current_factor_for_eppi0_bkg);
        } else {
            std::cout << "[current_dependence] Legacy post-binning current correction disabled; "
                      << "total_counts.cpp will fill corrected DATA/MC yields event-by-event from "
                      << options.response_model_json << std::endl;
        }

        write_csv_atomic(csv_path, csv);

        std::cout << "[current_dependence] Updated current-efficiency calibration columns and "
                  << "unity eppi0 normalization fallback columns in: " << csv_path << std::endl;
        if (options.apply_legacy_binned_current_corrections) {
            std::cout << "[current_dependence] Legacy mode also filled normalized DATA and "
                      << "current-corrected MC yield columns." << std::endl;
        } else {
            std::cout << "[current_dependence] Nominal corrected yield columns are deferred to "
                      << "total_counts.cpp event-level accumulation." << std::endl;
        }

        return true;
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return false;
    }
}