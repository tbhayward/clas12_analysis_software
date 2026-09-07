// cross_sections.cpp
// -----------------------------------------------------------------------------
// Cross section computation and plotting for DVCS ep -> ep gamma.
//
// Refactored workflow contract:
//
//   - All event-level yield corrections are applied upstream:
//       * current-efficiency correction,
//       * eppi0/AAOGEN cross-section normalization,
//       * pi0 contamination subtraction,
//       * acceptance unfolding.
//
//   - This module does NOT read imports/efficiency.json and does NOT apply any
//     additional efficiency or normalization scale.
//
//   - This module reads the already-unfolded yields from:
//       acceptance corrected yield, ep->epg, exp, <label>, <helicity>
//
//   - For single run-period labels, the cross section is computed directly from
//     that period's acceptance-corrected yield and that period's luminosity.
//
//   - [FIX] For combined labels, this module now recomputes the combined
//     acceptance-corrected yield row-by-row from the member-period
//     acceptance-corrected yields, using the same member-period validity mask
//     for both the numerator and the luminosity denominator. This avoids a
//     combined cross section being pulled by a signed/negative member-period
//     yield whose standalone cross section would not be reported.
//
//   - It reads correction factors already stored in the CSV:
//       Frad, <energy>
//       Fbin, <energy>
//       bin_volume, <energy>
//
//   - It computes the physical four-fold cross section in nb/(GeV^4 deg):
//       L_int[nb^-1] = Q[nC] * 1.316875 nb^-1/nC
//       sigma = Y_unfolded * Frad * Fbin / (L_int * bin_volume)
//
//   - It writes:
//       acceptance corrected yield, ep->epg, exp, <combined label>, <helicity>
//       cross sections, ep->epg, exp, <label>, <helicity>
//
//   - It also fills integrated luminosity columns in the CSV. For combined
//     labels, luminosity remains row-dependent and includes only member periods
//     that pass the same validity mask used in the combined yield numerator.
//
// Plotting and theory JSON generation are preserved from the previous version.
// -----------------------------------------------------------------------------

#include "cross_sections.h"
#include "model_predictions.h"

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
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <nlohmann/json.hpp>

// ROOT includes
#include <TCanvas.h>
#include <TError.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TPad.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TString.h>
#include <TH1.h>

namespace fs = std::filesystem;
using json   = nlohmann::json;

using Range = std::pair<double, double>;

// -----------------------------------------------------------------------------
// Configuration
// -----------------------------------------------------------------------------

// Frad/Fbin are read directly from Lee's imports/all_bin_v3.csv. The same
// Lee Frad/Fbin values are written into both the 10.6 GeV and 10.2 GeV
// pass-2 CSV columns before the cross sections are computed. The bin volume
// used in the denominator is the phase-space-allowed value already computed
// by bin_volume.cpp and stored in the pass-2 CSV.

static std::string canonical_period_dir(const std::string &label) {
    if (label == "Fa18 Inb")      return "Fa18_Inb";
    if (label == "Fa18 Out")      return "Fa18_Out";
    if (label == "Fa18 Inb Supp") return "Fa18_Inb_Supp";
    if (label == "Sp18 Inb")      return "Sp18_Inb";
    if (label == "Sp18 Out")      return "Sp18_Out";
    if (label == "Sp19 Inb")      return "Sp19_Inb";
    if (label == "Fa18")          return "Fa18";
    if (label == "Sp18")          return "Sp18";
    if (label == "10.6 GeV")      return "10.6_GeV";
    if (label == "10.2 GeV")      return "10.2_GeV";

    std::string out = label;
    std::replace(out.begin(), out.end(), ' ', '_');
    return out;
}

static std::string yield_label_for(const std::string &label) {
    if (label == "10.6 GeV") return "2018 (10.6 GeV)";
    if (label == "10.2 GeV") return "Sp19 Inb";
    return label;
}

static double beam_energy_for_label(const std::string &label) {
    if (label == "Sp19 Inb" || label == "10.2 GeV") {
        return 10.2;
    }

    return 10.6;
}

// -----------------------------------------------------------------------------
// Theory/model helper wrappers
// -----------------------------------------------------------------------------

static Helicity helicity_from_string(const std::string &h) {
    if (h == "pos") return Helicity::Plus;
    if (h == "neg") return Helicity::Minus;
    return Helicity::Unpol;
}

static double eval_bh_xs(double Ebeam,
                         double xB,
                         double Q2,
                         double t_pos,
                         double phi_rad) {
    const double phi_deg = phi_rad * 180.0 / M_PI;
    return vgg_bh_only(xB, Q2, t_pos, phi_deg, Ebeam);
}

static double eval_km_xs(double Ebeam,
                         double xB,
                         double Q2,
                         double t_pos,
                         double phi_rad,
                         const std::string &h) {
    const double phi_deg = phi_rad * 180.0 / M_PI;
    const Helicity hel = helicity_from_string(h);
    return km15_xs(xB, Q2, t_pos, phi_deg, Ebeam, hel);
}

static double eval_vgg_xs(double Ebeam,
                          double xB,
                          double Q2,
                          double t_pos,
                          double phi_rad,
                          const std::string &h) {
    const double phi_deg = phi_rad * 180.0 / M_PI;
    const Helicity hel = helicity_from_string(h);
    return vgg_xs(xB, Q2, t_pos, phi_deg, Ebeam, hel);
}

static void ensure_dir(const fs::path &p) {
    fs::create_directories(p);
}

// -----------------------------------------------------------------------------
// Basic CSV helpers
// -----------------------------------------------------------------------------

static std::vector<std::string> split_csv_line(const std::string &line) {
    std::vector<std::string> out;
    std::string field;
    bool in_quotes = false;

    for (size_t i = 0; i < line.size(); ++i) {
        char c = line[i];

        if (c == '"') {
            in_quotes = !in_quotes;
            field.push_back(c);
        } else if (c == ',' && !in_quotes) {
            out.push_back(field);
            field.clear();
        } else {
            field.push_back(c);
        }
    }

    out.push_back(field);
    return out;
}

static std::string trim(const std::string &s) {
    size_t b = 0;
    while (b < s.size() && std::isspace((unsigned char)s[b])) ++b;

    size_t e = s.size();
    while (e > b && std::isspace((unsigned char)s[e - 1])) --e;

    return s.substr(b, e - b);
}

static std::string unquote(const std::string &s) {
    if (s.size() >= 2 && s.front() == '"' && s.back() == '"') {
        std::string inner = s.substr(1, s.size() - 2);
        std::string out;

        for (size_t i = 0; i < inner.size(); ++i) {
            if (inner[i] == '"' && i + 1 < inner.size() && inner[i + 1] == '"') {
                out.push_back('"');
                ++i;
            } else {
                out.push_back(inner[i]);
            }
        }

        return out;
    }

    return s;
}

static std::string quote_if_needed(const std::string &s) {
    bool need = false;

    for (char c : s) {
        if (c == ',' || c == '"' || std::isspace((unsigned char)c)) {
            need = true;
            break;
        }
    }

    if (!need) return s;

    std::string out = "\"";

    for (char c : s) {
        if (c == '"') out += "\"\"";
        else out += c;
    }

    out += "\"";
    return out;
}

static std::string join_csv_line(const std::vector<std::string> &fields) {
    std::ostringstream oss;

    for (size_t i = 0; i < fields.size(); ++i) {
        if (i > 0) oss << ",";
        oss << quote_if_needed(fields[i]);
    }

    return oss.str();
}

static int find_col(const std::vector<std::string> &header,
                    const std::string &target) {
    for (size_t i = 0; i < header.size(); ++i) {
        if (trim(unquote(header[i])) == target) {
            return (int)i;
        }
    }

    throw std::runtime_error("Missing required column: \"" + target + "\"");
}

static int find_col_optional(const std::vector<std::string> &header,
                             const std::string &target) {
    for (size_t i = 0; i < header.size(); ++i) {
        if (trim(unquote(header[i])) == target) {
            return (int)i;
        }
    }

    return -1;
}

// -----------------------------------------------------------------------------
// Triple helpers (value, stat, sys) from cross_sections.h
// -----------------------------------------------------------------------------

// Helper: strip all outer quote layers (CSV and nested quotes)
static std::string strip_all_outer_quotes(std::string s) {
    s = unquote(s);
    s = trim(s);

    bool changed = true;

    while (changed && s.size() >= 2) {
        changed = false;

        char first = s.front();
        char last  = s.back();

        if ((first == '"' && last == '"') ||
            (first == '\'' && last == '\'')) {
            s = s.substr(1, s.size() - 2);
            s = trim(s);
            changed = true;
        }
    }

    return s;
}

static Triple parse_tuple3(const std::string &cell) {
    Triple out{0.0, 0.0, 0.0};

    // Keep the original quote/parenthesis tolerance, but parse the fixed
    // three-number payload directly. This avoids constructing a token vector
    // and three temporary substrings for every CSV tuple.
    std::string s = strip_all_outer_quotes(cell);
    s = trim(s);

    if (s.empty()) return out;

    if (s.front() == '(' && s.back() == ')') {
        s = s.substr(1, s.size() - 2);
        s = trim(s);
        if (s.empty()) return out;
    }

    const char *cursor = s.c_str();
    char *end = nullptr;

    auto parse_next = [&](double &value) -> bool {
        while (*cursor != '\0' && std::isspace(static_cast<unsigned char>(*cursor))) {
            ++cursor;
        }

        if (*cursor == '\0') return false;

        value = std::strtod(cursor, &end);
        if (end == cursor) return false;
        cursor = end;

        while (*cursor != '\0' && std::isspace(static_cast<unsigned char>(*cursor))) {
            ++cursor;
        }

        if (*cursor == ',') ++cursor;
        return true;
    };

    (void)parse_next(out.value);
    (void)parse_next(out.stat);
    (void)parse_next(out.sys);
    return out;
}

static std::string tuple3_to_cell(double value, double stat, double sys) {
    std::ostringstream oss;

    oss << "("
        << std::setprecision(10) << value << ", "
        << std::setprecision(10) << stat  << ", "
        << std::setprecision(10) << sys   << ")";

    return oss.str();
}

static Triple add_triples_quadrature_errors(const std::vector<Triple> &terms) {
    Triple out{0.0, 0.0, 0.0};

    double stat2 = 0.0;
    double sys2  = 0.0;

    for (const auto &t : terms) {
        out.value += t.value;

        if (std::isfinite(t.stat)) {
            stat2 += t.stat * t.stat;
        }

        if (std::isfinite(t.sys)) {
            sys2 += t.sys * t.sys;
        }
    }

    out.stat = std::sqrt(stat2);
    out.sys  = std::sqrt(sys2);

    return out;
}

// -----------------------------------------------------------------------------
// Luminosity helpers
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// RGA charge / luminosity inputs
// -----------------------------------------------------------------------------
//
// The single authoritative charge source for the production analysis is
//
//     imports/integrated_luminosity/global.csv
//
// Column convention:
//   column 1 = run number
//   column 2 = QADB-filtered RUN::Scaler accumulated charge [nC]
//   column 3 = positive-helicity accumulated charge [nC]
//   column 4 = negative-helicity accumulated charge [nC]
//
// The exact final Pass-2 run selections below are the same selections used by
// the run-selection/current-dependence analysis and documented in the Pass-2
// analysis note. Low-current Sp19 runs 6616 and 6618 are not included.
//
// LumiMap keeps the historical in-memory convention:
//   value = selected unpolarized accumulated charge [nC]
//   stat  = selected positive-helicity accumulated charge [nC]
//   sys   = selected negative-helicity accumulated charge [nC]
//
// Physical integrated luminosity is formed when the cross section is evaluated.
// For the 5-cm liquid-hydrogen target:
//
//   L_int = Q[nC] * 1.316875 nb^{-1}/nC
//
// equivalent to 1316.875 pb^{-1}/mC.
// -----------------------------------------------------------------------------

static constexpr double RGA_LUMINOSITY_NB_INV_PER_NC = 1.316875;
static constexpr double RGA_LUMINOSITY_PB_INV_PER_MC = 1316.875;

struct GlobalChargeRow {
    double total_nC = 0.0;
    double pos_nC = 0.0;
    double neg_nC = 0.0;
};

static const std::map<std::string, std::vector<int>>& final_pass2_run_selection() {
    static const std::map<std::string, std::vector<int>> runs = {
        {"Sp18 Inb", {
            3306, 3307, 3315, 3333, 3353, 3359, 3361, 3363, 3378, 3379, 3384, 3389,
            3390, 3403, 3405, 3406, 3407, 3409, 3411, 3421, 3422, 3429, 3431, 3432,
            3433, 3434, 3435, 3436, 3441, 3442, 3459, 3460, 3461, 3462, 3463, 3464,
            3465, 3466, 3467, 3469, 3480, 3482, 3484, 3485, 3488, 3492, 3493, 3501,
            3506, 3507, 3512, 3513, 3517, 3518, 3519, 3520, 3521, 3522, 3542, 3699,
            3700, 3702, 3705, 3708, 3711, 3719, 3720, 3722, 3738, 3739, 3741, 3748,
            3750, 3752, 3771, 3789, 3790, 3791, 3795, 3796, 3797, 3798, 3799, 3803,
            3804, 3806, 3816, 4003, 4013, 4014, 4015, 4016, 4017, 4021, 4022, 4025,
            4026, 4028, 4030, 4032, 4033, 4037, 4038, 4039, 4041, 4044, 4045, 4050,
            4053, 4054, 4055, 4058, 4060, 4061, 4067, 4068, 4069, 4070, 4071, 4073,
            4075, 4078, 4080, 4081, 4082, 4083, 4085, 4089, 4090, 4091, 4092, 4093,
            4094, 4095, 4096, 4097, 4099, 4100, 4103, 4104, 4110, 4112, 4113, 4114,
            4115, 4139, 4143, 4144, 4147, 4148, 4151, 4152, 4154, 4155, 4156, 4157,
            4158, 4161, 4164, 4165, 4166, 4167, 4169, 4174, 4180, 4181, 4182, 4184,
            4189, 4190, 4192, 4193, 4194, 4201, 4202, 4203, 4204, 4206, 4208, 4211,
            4213, 4218, 4219, 4220, 4221, 4223, 4224, 4226, 4229, 4243, 4244, 4245,
            4247, 4248, 4250, 4253, 4254, 4256, 4262, 4263, 4264, 4309, 4311, 4312,
            4313, 4314, 4315, 4316, 4320, 4321, 4322, 4323, 4324
        }},
        {"Sp18 Out", {
            3261, 3262, 3266, 3269, 3270, 3282, 3288, 3874, 3875, 3878, 3880, 3881,
            3883, 3884, 3885, 3888, 3889, 3891, 3893, 3898, 3903, 3905, 3907, 3908,
            3910, 3911, 3912, 3913, 3915, 3916, 3917, 3919, 3920, 3921, 3924, 3926,
            3928, 3930, 3932, 3933, 3934, 3936, 3938, 3939, 3940, 3941, 3943, 3944,
            3945, 3946, 3948, 3949, 3950, 3954, 3959, 3963, 3964, 3969, 3970, 3973,
            3975, 3982, 3985, 3986, 3987
        }},
        {"Fa18 Inb", {
            5335, 5339, 5340, 5341, 5342, 5343, 5344, 5032, 5036, 5038, 5039, 5040,
            5041, 5043, 5045, 5046, 5047, 5051, 5052, 5053, 5116, 5117, 5119, 5120,
            5124, 5125, 5126, 5127, 5128, 5129, 5130, 5139, 5153, 5158, 5159, 5160,
            5162, 5163, 5164, 5165, 5166, 5167, 5168, 5169, 5180, 5181, 5182, 5183,
            5190, 5191, 5193, 5195, 5196, 5197, 5198, 5199, 5200, 5201, 5202, 5203,
            5204, 5205, 5206, 5208, 5211, 5212, 5215, 5216, 5219, 5220, 5221, 5222,
            5223, 5230, 5231, 5232, 5233, 5234, 5235, 5237, 5238, 5239, 5248, 5249,
            5252, 5253, 5257, 5258, 5259, 5261, 5262, 5303, 5304, 5305, 5306, 5307,
            5310, 5311, 5315, 5317, 5318, 5319, 5320, 5323, 5324, 5333, 5334, 5336,
            5346, 5347, 5349, 5351, 5354, 5355, 5367, 5356, 5357, 5358, 5359, 5360,
            5361, 5362, 5366, 5368, 5369, 5372, 5373, 5374, 5375, 5376, 5377, 5378,
            5379, 5380, 5381, 5382, 5383, 5386, 5390, 5391, 5392, 5393, 5398, 5400,
            5401, 5403, 5404, 5406, 5407
        }},
        {"Fa18 Out", {
            5444, 5423, 5424, 5425, 5426, 5428, 5429, 5430, 5432, 5434, 5435, 5436,
            5437, 5438, 5440, 5441, 5442, 5445, 5447, 5448, 5449, 5450, 5451, 5452,
            5453, 5454, 5455, 5460, 5464, 5465, 5466, 5467, 5468, 5469, 5470, 5471,
            5472, 5473, 5474, 5475, 5476, 5478, 5479, 5480, 5481, 5482, 5483, 5485,
            5486, 5487, 5495, 5496, 5497, 5498, 5499, 5500, 5504, 5505, 5507, 5516,
            5517, 5518, 5519, 5520, 5521, 5522, 5523, 5524, 5525, 5526, 5527, 5528,
            5530, 5532, 5533, 5534, 5535, 5536, 5537, 5538, 5540, 5541, 5543, 5544,
            5545, 5546, 5547, 5548, 5549, 5550, 5551, 5552, 5555, 5556, 5557, 5558,
            5559, 5562, 5567, 5569, 5570, 5571, 5572, 5573, 5574, 5577, 5578, 5591,
            5592, 5594, 5597, 5598, 5600, 5601, 5602, 5603, 5604, 5606, 5607, 5611,
            5612, 5613, 5614, 5615, 5616, 5617, 5618, 5619, 5621, 5623, 5624, 5625,
            5626, 5627, 5628, 5629, 5630, 5631, 5632, 5633, 5635, 5637, 5638, 5639,
            5641, 5643, 5644, 5645, 5646, 5647, 5648, 5649, 5650, 5651, 5652, 5654,
            5655, 5656, 5662, 5663, 5664, 5665, 5666
        }},
        {"Sp19 Inb", {
            6619, 6620, 6636, 6637, 6638, 6639, 6640, 6642, 6645, 6647, 6648, 6650,
            6651, 6652, 6654, 6655, 6656, 6657, 6658, 6660, 6661, 6662, 6663, 6664,
            6665, 6666, 6667, 6668, 6669, 6670, 6672, 6673, 6675, 6676, 6677, 6678,
            6680, 6682, 6683, 6684, 6685, 6687, 6688, 6689, 6691, 6692, 6693, 6694,
            6695, 6696, 6697, 6698, 6699, 6704, 6705, 6706, 6707, 6708, 6709, 6710,
            6711, 6712, 6713, 6714, 6715, 6716, 6717, 6718, 6719, 6729, 6730, 6731,
            6732, 6733, 6736, 6737, 6738, 6739, 6740, 6741, 6742, 6743, 6744, 6746,
            6747, 6748, 6749, 6750, 6753, 6754, 6755, 6756, 6757, 6759, 6760, 6762,
            6763, 6764, 6765, 6767, 6768, 6769, 6779, 6780, 6781, 6783
        }}
    };
    return runs;
}

static std::map<int, GlobalChargeRow> load_global_charge_csv(const std::string &path) {
    std::ifstream ifs(path);
    if (!ifs) {
        throw std::runtime_error(
            "[cross_sections] cannot open authoritative charge CSV: " + path
        );
    }

    std::map<int, GlobalChargeRow> out;
    std::string line;
    int line_number = 0;

    while (std::getline(ifs, line)) {
        ++line_number;
        const std::string s = trim(line);
        if (s.empty() || s[0] == '#') continue;

        const std::vector<std::string> fields = split_csv_line(s);
        if (fields.size() < 4) {
            std::ostringstream ss;
            ss << "[cross_sections] malformed charge row at " << path
               << ":" << line_number << " (need at least 4 columns)";
            throw std::runtime_error(ss.str());
        }

        int run = 0;
        double total = 0.0;
        double pos = 0.0;
        double neg = 0.0;

        try {
            run   = std::stoi(trim(unquote(fields[0])));
            total = std::stod(trim(unquote(fields[1])));
            pos   = std::stod(trim(unquote(fields[2])));
            neg   = std::stod(trim(unquote(fields[3])));
        } catch (const std::exception &) {
            std::ostringstream ss;
            ss << "[cross_sections] invalid numeric value at " << path
               << ":" << line_number;
            throw std::runtime_error(ss.str());
        }

        if (out.count(run)) {
            std::ostringstream ss;
            ss << "[cross_sections] duplicate run " << run << " in " << path;
            throw std::runtime_error(ss.str());
        }

        out[run] = GlobalChargeRow{total, pos, neg};
    }

    if (out.empty()) {
        throw std::runtime_error(
            "[cross_sections] no charge records loaded from " + path
        );
    }

    return out;
}

static Triple sum_selected_charge_for_period(
    const std::string &period,
    const std::map<int, GlobalChargeRow> &charge_rows) {

    const auto &selection = final_pass2_run_selection();
    const auto it_sel = selection.find(period);

    if (it_sel == selection.end()) {
        throw std::runtime_error(
            "[cross_sections] no final Pass-2 run selection defined for " + period
        );
    }

    Triple out{0.0, 0.0, 0.0};

    for (const int run : it_sel->second) {
        const auto it = charge_rows.find(run);
        if (it == charge_rows.end()) {
            std::ostringstream ss;
            ss << "[cross_sections] selected run " << run
               << " (" << period << ") is missing from global.csv";
            throw std::runtime_error(ss.str());
        }

        out.value += it->second.total_nC;
        out.stat  += it->second.pos_nC;
        out.sys   += it->second.neg_nC;
    }

    return out;
}

static double integrated_luminosity_nb_inv(double charge_nC) {
    return charge_nC * RGA_LUMINOSITY_NB_INV_PER_NC;
}

static double integrated_luminosity_pb_inv(double charge_nC) {
    return charge_nC * 1.0e-6 * RGA_LUMINOSITY_PB_INV_PER_MC;
}

LumiMap build_lumi_map() {
    LumiBuildOptions options;
    return build_lumi_map(options);
}

LumiMap build_lumi_map(const LumiBuildOptions &options) {
    const std::string charge_csv = options.charge_csv_path;

    std::cout << "[cross_sections] Charge source: " << charge_csv << "\n"
              << "[cross_sections] Unpolarized normalization: column 2 "
              << "(QADB-filtered RUN::Scaler charge) over the final Pass-2 run selection.\n";

    LumiMap m;

    try {
        const std::map<int, GlobalChargeRow> rows =
            load_global_charge_csv(charge_csv);

        for (const auto &period : std::vector<std::string>{
                 "Fa18 Inb", "Fa18 Out", "Sp18 Inb", "Sp18 Out", "Sp19 Inb"}) {
            m[period] = sum_selected_charge_for_period(period, rows);
        }
    } catch (const std::exception &e) {
        std::cerr << "[cross_sections] FATAL in build_lumi_map: "
                  << e.what() << "\n";
        throw;
    }

    m["Fa18 Inb Supp"] = Triple{0.0, 0.0, 0.0};

    auto sum_labels = [&](const std::vector<std::string> &labels) {
        Triple out{0.0, 0.0, 0.0};
        for (const auto &label : labels) {
            auto it = m.find(label);
            if (it == m.end()) continue;
            out.value += it->second.value;
            out.stat  += it->second.stat;
            out.sys   += it->second.sys;
        }
        return out;
    };

    m["Fa18"]     = sum_labels({"Fa18 Inb", "Fa18 Out"});
    m["Sp18"]     = sum_labels({"Sp18 Inb", "Sp18 Out"});
    m["10.6 GeV"] = sum_labels({"Fa18 Inb", "Fa18 Out", "Sp18 Inb", "Sp18 Out"});
    m["10.2 GeV"] = sum_labels({"Sp19 Inb"});

    std::cout << std::fixed << std::setprecision(6)
              << "[cross_sections] Selected accumulated charge from global.csv:\n";

    for (const auto &period : std::vector<std::string>{
             "Sp18 Inb", "Sp18 Out", "Fa18 Inb", "Fa18 Out", "Sp19 Inb"}) {
        const Triple &q = m.at(period);
        std::cout << "  " << std::setw(9) << std::left << period
                  << "  Q=" << std::setw(10) << std::right
                  << q.value / 1.0e6 << " mC"
                  << "  L_int=" << integrated_luminosity_pb_inv(q.value) / 1000.0
                  << " fb^-1\n";
    }

    std::cout << std::defaultfloat << std::setprecision(6);
    return m;
}

// -----------------------------------------------------------------------------
// Theory JSON generation (xs_phi_all.json)
// -----------------------------------------------------------------------------

static bool write_theory_json_for_energy(const std::string &csv_main,
                                         const std::string &theory_json_root,
                                         double Ebeam,
                                         const std::string &energy_label) {
    std::vector<double> phi_deg;
    phi_deg.reserve(38);

    const int    N       = 38;
    const double phi_min = 0.1;
    const double phi_max = 359.9;
    const double step    = (phi_max - phi_min) / (double)(N - 1);

    for (int i = 0; i < N; ++i) {
        double phi = phi_min + step * (double)i;

        if (i == N - 1) {
            phi = phi_max;
        }

        phi_deg.push_back(phi);
    }

    const int n_phi = (int)phi_deg.size();

    std::ifstream ifs(csv_main);

    if (!ifs) {
        std::cerr << "[cross_sections] FATAL: cannot open " << csv_main
                  << " for theory JSON generation (energy " << energy_label
                  << ").\n";
        return false;
    }

    std::vector<std::string> lines;
    std::string line;

    while (std::getline(ifs, line)) {
        lines.push_back(line);
    }

    ifs.close();

    if (lines.empty()) {
        std::cerr << "[cross_sections] FATAL: CSV " << csv_main
                  << " is empty in write_theory_json_for_energy.\n";
        return false;
    }

    std::vector<std::string> header = split_csv_line(lines[0]);

    int c_xb_min = -1;
    int c_xb_max = -1;
    int c_q2_min = -1;
    int c_q2_max = -1;
    int c_t_min  = -1;
    int c_t_max  = -1;

    try {
        c_xb_min = find_col(header, "xBmin");
        c_xb_max = find_col(header, "xBmax");
        c_q2_min = find_col(header, "Q2min");
        c_q2_max = find_col(header, "Q2max");
        c_t_min  = find_col(header, "t_abs_min");
        c_t_max  = find_col(header, "t_abs_max");
    } catch (const std::exception &e) {
        std::cerr << "[cross_sections] FATAL: " << e.what()
                  << " in write_theory_json_for_energy.\n";
        return false;
    }

    json j;
    j["phi_deg"] = phi_deg;
    json rows_json = json::object();

    const size_t n_rows = lines.size();

    std::cout << "[cross_sections] Generating theory JSON for energy \""
              << energy_label << "\" (Ebeam=" << Ebeam
              << ") for " << (n_rows > 0 ? n_rows - 1 : 0)
              << " data rows.\n";

    int next_pct = 1;

    for (size_t row = 1; row < n_rows; ++row) {
        if (lines[row].empty()) continue;

        if (n_rows > 1 && next_pct <= 100) {
            double frac = 100.0 * (double)row / (double)(n_rows - 1);

            if (frac >= next_pct) {
                std::cout << "[cross_sections] theory JSON ("
                          << energy_label << "): ~"
                          << next_pct << "% of rows processed (row "
                          << row << " / " << (n_rows - 1) << ")\n";
                next_pct += 10;
            }
        }

        std::vector<std::string> fields = split_csv_line(lines[row]);

        if (fields.size() != header.size()) {
            std::cerr << "[cross_sections] WARNING: row " << row
                      << " has wrong number of fields, skipping.\n";
            continue;
        }

        double xbmin = std::atof(trim(unquote(fields[c_xb_min])).c_str());
        double xbmax = std::atof(trim(unquote(fields[c_xb_max])).c_str());
        double q2min = std::atof(trim(unquote(fields[c_q2_min])).c_str());
        double q2max = std::atof(trim(unquote(fields[c_q2_max])).c_str());
        double tmin  = std::atof(trim(unquote(fields[c_t_min])).c_str());
        double tmax  = std::atof(trim(unquote(fields[c_t_max])).c_str());

        double xB_mid = 0.5 * (xbmin + xbmax);
        double Q2_mid = 0.5 * (q2min + q2max);
        double t_mid  = 0.5 * (tmin  + tmax);

        if (!(xB_mid > 0.0) || !(Q2_mid > 0.0) || !(t_mid > 0.0)) {
            continue;
        }

        std::vector<double> bh_unpol(n_phi);
        std::vector<double> bh_pos(n_phi);
        std::vector<double> bh_neg(n_phi);

        std::vector<double> km_unpol(n_phi);
        std::vector<double> km_pos(n_phi);
        std::vector<double> km_neg(n_phi);

        std::vector<double> vgg_unpol(n_phi);
        std::vector<double> vgg_pos(n_phi);
        std::vector<double> vgg_neg(n_phi);

        for (int i = 0; i < n_phi; ++i) {
            double phideg = phi_deg[i];
            double phirad = phideg * M_PI / 180.0;

            double bh = eval_bh_xs(Ebeam, xB_mid, Q2_mid, t_mid, phirad);
            bh_unpol[i] = bh;
            bh_pos[i]   = bh;
            bh_neg[i]   = bh;

            km_unpol[i] = eval_km_xs(Ebeam, xB_mid, Q2_mid, t_mid, phirad, "unpol");
            km_pos[i]   = eval_km_xs(Ebeam, xB_mid, Q2_mid, t_mid, phirad, "pos");
            km_neg[i]   = eval_km_xs(Ebeam, xB_mid, Q2_mid, t_mid, phirad, "neg");

            vgg_unpol[i] = eval_vgg_xs(Ebeam, xB_mid, Q2_mid, t_mid, phirad, "unpol");
            vgg_pos[i]   = eval_vgg_xs(Ebeam, xB_mid, Q2_mid, t_mid, phirad, "pos");
            vgg_neg[i]   = eval_vgg_xs(Ebeam, xB_mid, Q2_mid, t_mid, phirad, "neg");
        }

        json bh_json;
        json km_json;
        json vgg_json;

        bh_json["unpol"] = bh_unpol;
        bh_json["pos"]   = bh_pos;
        bh_json["neg"]   = bh_neg;

        km_json["unpol"] = km_unpol;
        km_json["pos"]   = km_pos;
        km_json["neg"]   = km_neg;

        vgg_json["unpol"] = vgg_unpol;
        vgg_json["pos"]   = vgg_pos;
        vgg_json["neg"]   = vgg_neg;

        json row_json;
        row_json["BH"]  = bh_json;
        row_json["KM"]  = km_json;
        row_json["VGG"] = vgg_json;

        rows_json[std::to_string(row)] = std::move(row_json);
    }

    j["rows"] = std::move(rows_json);

    fs::path dir  = fs::path(theory_json_root) / canonical_period_dir(energy_label);
    ensure_dir(dir);

    fs::path file = dir / "xs_phi_all.json";

    std::ofstream ofs(file);

    if (!ofs) {
        std::cerr << "[cross_sections] FATAL: cannot open "
                  << file.string()
                  << " for writing theory JSON.\n";
        return false;
    }

    ofs << std::setw(2) << j << "\n";
    ofs.close();

    std::cout << "[cross_sections] Wrote theory JSON xs_phi_all.json for energy \""
              << energy_label << "\" at " << file.string() << "\n";

    return true;
}

bool regenerate_theory_jsons(const std::string &csv_main,
                             const std::string &theory_json_root) {
    bool ok_106 = write_theory_json_for_energy(csv_main, theory_json_root,
                                               10.6, "10.6 GeV");

    bool ok_102 = write_theory_json_for_energy(csv_main, theory_json_root,
                                               10.2, "10.2 GeV");

    if (!ok_106 || !ok_102) {
        std::cerr << "[cross_sections] ERROR: regenerate_theory_jsons failed for "
                  << "10.6 (ok=" << ok_106 << ") or 10.2 (ok=" << ok_102 << ").\n";
        return false;
    }

    std::cout << "[cross_sections] regenerate_theory_jsons completed.\n";
    return true;
}

// -----------------------------------------------------------------------------
// Cross section computation
// -----------------------------------------------------------------------------

static bool is_combined_label(const std::string &L) {
    return (L == "Fa18" || L == "Sp18" || L == "10.6 GeV");
}

static std::vector<std::string> combined_members_for(const std::string &L) {
    if (L == "Fa18") return {"Fa18 Inb", "Fa18 Out"};
    if (L == "Sp18") return {"Sp18 Inb", "Sp18 Out"};
    if (L == "10.6 GeV") return {"Fa18 Inb", "Fa18 Out", "Sp18 Inb", "Sp18 Out"};

    return {};
}

static const std::vector<std::string> kAllHelicities = {"unpol", "pos", "neg"};
static const std::vector<std::string> kUnpolarizedOnlyHelicities = {"unpol"};

static bool has_helicity_resolved_cross_sections(const std::string &label) {
    if (label == "Sp18 Inb" || label == "Sp18 Out" ||
        label == "Sp18" || label == "10.6 GeV") {
        return false;
    }

    return true;
}

static const std::vector<std::string>& cross_section_helicities_for_label(const std::string &label) {
    if (has_helicity_resolved_cross_sections(label)) {
        return kAllHelicities;
    }

    return kUnpolarizedOnlyHelicities;
}

static std::string lumi_col_for_label(const std::string &L) {
    if (L == "Fa18 Inb") return "integrated luminosity, Fa18 Inb (nC)";
    if (L == "Fa18 Out") return "integrated luminosity, Fa18 Out (nC)";
    if (L == "Sp19 Inb") return "integrated luminosity, Sp19 Inb (nC)";
    if (L == "Sp18 Inb") return "integrated luminosity, Sp18 Inb (nC)";
    if (L == "Sp18 Out") return "integrated luminosity, Sp18 Out (nC)";
    if (L == "Fa18") return "integrated luminosity, Fa18 (nC)";
    if (L == "Sp18") return "integrated luminosity, Sp18 (nC)";
    if (L == "10.6 GeV") return "integrated luminosity, 10.6 GeV (nC)";

    return "";
}

static double lumi_component_for_helicity(const Triple &Lumi,
                                          const std::string &helicity) {
    if (helicity == "unpol") return Lumi.value;
    if (helicity == "pos")   return Lumi.stat;
    if (helicity == "neg")   return Lumi.sys;

    return 0.0;
}

static double ratio_rel2(double value, double err) {
    if (value == 0.0 || !std::isfinite(value) || !std::isfinite(err) || err <= 0.0) {
        return 0.0;
    }

    const double r = err / value;
    return r * r;
}

// -----------------------------------------------------------------------------
// Lee correction-factor import
// -----------------------------------------------------------------------------

struct LeeCorrections {
    Triple frad{0.0, 0.0, 0.0};
    Triple fbin{0.0, 0.0, 0.0};
};

static std::string first_nonempty_bin_index(const std::vector<std::string> &fields,
                                            int c_bin_index,
                                            int one_based_row_number) {
    if (c_bin_index >= 0 && c_bin_index < (int)fields.size()) {
        const std::string v = trim(unquote(fields[c_bin_index]));

        if (!v.empty()) return v;
    }

    std::ostringstream ss;
    ss << one_based_row_number;
    return ss.str();
}

static double parse_required_number(const std::string &cell,
                                    const std::string &label,
                                    const std::string &source) {
    const std::string s = trim(unquote(cell));

    if (s.empty()) {
        std::ostringstream ss;
        ss << "[cross_sections] FATAL: empty numeric cell for " << label
           << " in " << source;
        throw std::runtime_error(ss.str());
    }

    char *endp = nullptr;
    const double v = std::strtod(s.c_str(), &endp);

    if (endp == s.c_str() || !std::isfinite(v)) {
        std::ostringstream ss;
        ss << "[cross_sections] FATAL: could not parse numeric cell for " << label
           << " from value '" << s << "' in " << source;
        throw std::runtime_error(ss.str());
    }

    return v;
}

static Triple parse_required_triple_cell(const std::vector<std::string> &fields,
                                         int col,
                                         const std::string &col_name,
                                         const std::string &context) {
    if (col < 0 || col >= (int)fields.size()) {
        std::ostringstream ss;
        ss << "[cross_sections] FATAL: column index out of range for "
           << col_name << " while reading " << context;
        throw std::runtime_error(ss.str());
    }

    const std::string s = trim(strip_all_outer_quotes(fields[col]));

    if (s.empty()) {
        std::ostringstream ss;
        ss << "[cross_sections] FATAL: empty required tuple cell for "
           << col_name << " while reading " << context;
        throw std::runtime_error(ss.str());
    }

    Triple out = parse_tuple3(fields[col]);

    if (!std::isfinite(out.value) || out.value <= 0.0) {
        std::ostringstream ss;
        ss << "[cross_sections] FATAL: required tuple cell for "
           << col_name << " is missing, malformed, or non-positive while reading "
           << context << ". Cell value was '" << s << "'.";
        throw std::runtime_error(ss.str());
    }

    if (!std::isfinite(out.stat)) out.stat = 0.0;
    if (!std::isfinite(out.sys))  out.sys  = 0.0;

    return out;
}

static std::map<std::string, LeeCorrections>
load_lee_corrections_by_bin_index(const std::string &lee_csv_path) {
    std::ifstream ifs(lee_csv_path);

    if (!ifs) {
        std::ostringstream ss;
        ss << "[cross_sections] FATAL: cannot open Lee correction CSV: "
           << lee_csv_path;
        throw std::runtime_error(ss.str());
    }

    std::string line;

    if (!std::getline(ifs, line)) {
        std::ostringstream ss;
        ss << "[cross_sections] FATAL: Lee correction CSV is empty: "
           << lee_csv_path;
        throw std::runtime_error(ss.str());
    }

    std::vector<std::string> header = split_csv_line(line);

    int c_bin_index = find_col_optional(header, "bin index");

    if (c_bin_index < 0) {
        c_bin_index = find_col_optional(header, "");
    }

    const int c_frad = find_col(header, "Frad");
    const int c_fbin = find_col(header, "Fbin");
    const int c_valid = find_col_optional(header, "valid bin");

    std::map<std::string, LeeCorrections> out;

    int input_row = 0;
    int kept_row = 0;

    while (std::getline(ifs, line)) {
        ++input_row;

        if (line.empty()) continue;

        std::vector<std::string> fields = split_csv_line(line);

        if (fields.size() < header.size()) {
            fields.resize(header.size());
        }

        if (fields.size() != header.size()) {
            std::ostringstream ss;
            ss << "[cross_sections] FATAL: Lee row width mismatch in "
               << lee_csv_path << " on input row " << input_row;
            throw std::runtime_error(ss.str());
        }

        if (c_valid >= 0) {
            const std::string valid_s = trim(unquote(fields[c_valid]));

            if (!(valid_s == "1" || valid_s == "1.0" ||
                  valid_s == "true" || valid_s == "TRUE")) {
                continue;
            }
        }

        const std::string bin_index = first_nonempty_bin_index(fields, c_bin_index, input_row);

        LeeCorrections c;
        c.frad.value = parse_required_number(fields[c_frad], "Frad", lee_csv_path);
        c.fbin.value = parse_required_number(fields[c_fbin], "Fbin", lee_csv_path);

        if (!(c.frad.value > 0.0) || !(c.fbin.value > 0.0)) {
            std::ostringstream ss;
            ss << "[cross_sections] FATAL: non-positive Lee correction for bin index "
               << bin_index << " in " << lee_csv_path
               << " (Frad=" << c.frad.value
               << ", Fbin=" << c.fbin.value << ")";
            throw std::runtime_error(ss.str());
        }

        if (!out.emplace(bin_index, c).second) {
            std::ostringstream ss;
            ss << "[cross_sections] FATAL: duplicate Lee bin index " << bin_index
               << " in " << lee_csv_path;
            throw std::runtime_error(ss.str());
        }

        ++kept_row;
    }

    std::cout << "[cross_sections] Loaded Lee Frad/Fbin correction factors from "
              << lee_csv_path << ": input rows=" << input_row
              << " valid rows=" << kept_row << "\n";

    return out;
}

static const LeeCorrections& lee_for_pass2_row(
    const std::map<std::string, LeeCorrections> &lee,
    const std::vector<std::string> &fields,
    int c_pass2_bin_index,
    size_t csv_row_number) {

    if (c_pass2_bin_index < 0 || c_pass2_bin_index >= (int)fields.size()) {
        throw std::runtime_error("[cross_sections] FATAL: pass-2 CSV is missing bin index column.");
    }

    const std::string bin_index = trim(unquote(fields[c_pass2_bin_index]));

    if (bin_index.empty()) {
        std::ostringstream ss;
        ss << "[cross_sections] FATAL: empty bin index on pass-2 CSV data row "
           << csv_row_number;
        throw std::runtime_error(ss.str());
    }

    auto it = lee.find(bin_index);

    if (it == lee.end()) {
        std::ostringstream ss;
        ss << "[cross_sections] FATAL: no Lee correction found for pass-2 bin index "
           << bin_index << " on CSV data row " << csv_row_number;
        throw std::runtime_error(ss.str());
    }

    return it->second;
}

bool compute_cross_sections(const std::string &csv_main,
                            const LumiMap &lumi_map) {
    return compute_cross_sections(csv_main, lumi_map, "imports/all_bin_v3.csv");
}

bool compute_cross_sections(const std::string &csv_main,
                            const LumiMap &lumi_map,
                            const std::string &lee_csv_path) {
    std::map<std::string, LeeCorrections> lee_corrections;

    try {
        lee_corrections = load_lee_corrections_by_bin_index(lee_csv_path);
    } catch (const std::exception &e) {
        std::cerr << e.what() << "\n";
        return false;
    }

    std::ifstream ifs(csv_main);

    if (!ifs) {
        std::cerr << "[cross_sections] ERROR: cannot open " << csv_main << " for reading.\n";
        return false;
    }

    std::vector<std::string> lines;
    std::string line;

    while (std::getline(ifs, line)) {
        lines.push_back(line);
    }

    ifs.close();

    if (lines.empty()) {
        std::cerr << "[cross_sections] ERROR: CSV " << csv_main << " is empty.\n";
        return false;
    }

    std::vector<std::string> header = split_csv_line(lines[0]);
    const size_t n_data_rows = lines.size() - 1;

    const int c_pass2_bin_index = find_col_optional(header, "bin index");

    if (c_pass2_bin_index < 0) {
        std::cerr << "[cross_sections] FATAL: missing pass-2 CSV column: bin index\n";
        return false;
    }

    const int c_vbin_106 = find_col_optional(header, "bin_volume, 10.6 GeV");
    const int c_vbin_102 = find_col_optional(header, "bin_volume, 10.2 GeV");
    const int c_cubic_vbin_106 = find_col_optional(header, "cubic bin_volume, 10.6 GeV");
    const int c_cubic_vbin_102 = find_col_optional(header, "cubic bin_volume, 10.2 GeV");
    const int c_frad_106 = find_col_optional(header, "Frad, 10.6 GeV");
    const int c_frad_102 = find_col_optional(header, "Frad, 10.2 GeV");
    const int c_fbin_106 = find_col_optional(header, "Fbin, 10.6 GeV");
    const int c_fbin_102 = find_col_optional(header, "Fbin, 10.2 GeV");

    if (c_vbin_106 < 0) {
        std::cerr << "[cross_sections] FATAL: missing phase-space bin volume column: bin_volume, 10.6 GeV\n";
        return false;
    }

    if (c_vbin_102 < 0) {
        std::cerr << "[cross_sections] FATAL: missing phase-space bin volume column: bin_volume, 10.2 GeV\n";
        return false;
    }

    if (c_cubic_vbin_106 < 0) {
        std::cerr << "[cross_sections] FATAL: missing diagnostic cubic bin volume column: cubic bin_volume, 10.6 GeV\n";
        return false;
    }

    if (c_cubic_vbin_102 < 0) {
        std::cerr << "[cross_sections] FATAL: missing diagnostic cubic bin volume column: cubic bin_volume, 10.2 GeV\n";
        return false;
    }

    if (c_frad_106 < 0) {
        std::cerr << "[cross_sections] FATAL: missing Frad column: Frad, 10.6 GeV\n";
        return false;
    }

    if (c_frad_102 < 0) {
        std::cerr << "[cross_sections] FATAL: missing Frad column: Frad, 10.2 GeV\n";
        return false;
    }

    if (c_fbin_106 < 0) {
        std::cerr << "[cross_sections] FATAL: missing Fbin column: Fbin, 10.6 GeV\n";
        return false;
    }

    if (c_fbin_102 < 0) {
        std::cerr << "[cross_sections] FATAL: missing Fbin column: Fbin, 10.2 GeV\n";
        return false;
    }

    const std::vector<std::string> base_periods = {
        "Fa18 Inb",
        "Fa18 Out",
        "Sp19 Inb",
        "Sp18 Inb",
        "Sp18 Out"
    };

    std::map<std::string, int> acc_col_idx;

    for (const auto &p : base_periods) {
        const std::string col = "acceptance, " + p;
        int idx = find_col_optional(header, col);

        if (idx < 0) {
            std::cerr << "[cross_sections] FATAL: missing required acceptance column: "
                      << col << "\n";
            return false;
        }

        acc_col_idx[p] = idx;
    }

    const std::vector<std::string> lumi_cols_required = {
        "integrated luminosity, Fa18 Inb (nC)",
        "integrated luminosity, Fa18 Out (nC)",
        "integrated luminosity, Sp19 Inb (nC)",
        "integrated luminosity, Sp18 Inb (nC)",
        "integrated luminosity, Sp18 Out (nC)",
        "integrated luminosity, Fa18 (nC)",
        "integrated luminosity, Sp18 (nC)",
        "integrated luminosity, 10.6 GeV (nC)"
    };

    std::map<std::string, int> lumi_col_idx;

    for (const auto &name : lumi_cols_required) {
        int idx = find_col_optional(header, name);

        if (idx < 0) {
            std::cerr << "[cross_sections] FATAL: missing luminosity column: "
                      << name << "\n";
            return false;
        }

        lumi_col_idx[name] = idx;
    }

    const std::vector<std::string> labels = {
        "Fa18 Inb",
        "Fa18 Out",
        "Sp18 Inb",
        "Sp18 Out",
        "Sp19 Inb",
        "Fa18",
        "Sp18",
        "10.6 GeV"
    };

    struct ColPair {
        int yield_idx = -1;
        int xs_idx = -1;
    };

    std::map<std::string, ColPair> colmap;

    for (const auto &L : labels) {
        const std::string YL = yield_label_for(L);

        for (const auto &h : cross_section_helicities_for_label(L)) {
            const std::string yield_col =
                "acceptance corrected yield, ep->epg, exp, " + YL + ", " + h;

            const std::string xs_col =
                "cross sections, ep->epg, exp, " + L + ", " + h;

            const int iy = find_col_optional(header, yield_col);
            const int ix = find_col_optional(header, xs_col);

            if (iy >= 0 && ix >= 0) {
                colmap[L + "|" + h] = ColPair{iy, ix};
            } else if (iy >= 0 && ix < 0) {
                std::cerr << "[cross_sections] FATAL: yield column exists but cross-section column is missing: "
                          << xs_col << "\n";
                return false;
            } else if (iy < 0 && ix >= 0) {
                std::cerr << "[cross_sections] FATAL: cross-section column exists but yield column is missing: "
                          << yield_col << "\n";
                return false;
            }
        }
    }

    // [FIX] Explicit map of single-period acceptance-corrected yield columns.
    // Combined labels are recomputed from these columns instead of trusting the
    // upstream pre-combined yield blindly.
    std::map<std::string, int> member_yield_col_idx;

    for (const auto &p : base_periods) {
        for (const auto &h : cross_section_helicities_for_label(p)) {
            const std::string yield_col =
                "acceptance corrected yield, ep->epg, exp, " + p + ", " + h;

            const int iy = find_col_optional(header, yield_col);

            if (iy >= 0) {
                member_yield_col_idx[p + "|" + h] = iy;
            }
        }
    }

    auto period_has_acceptance = [&](const std::string &period,
                                     const std::vector<std::string> &fields) -> bool {
        auto it = acc_col_idx.find(period);

        if (it == acc_col_idx.end()) return false;

        const Triple a = parse_tuple3(fields[it->second]);
        return a.value > 0.0 && std::isfinite(a.value);
    };

    auto member_period_lumi = [&](const std::string &period) -> Triple {
        auto it = lumi_map.find(period);

        if (it == lumi_map.end()) return Triple{0.0, 0.0, 0.0};

        return it->second;
    };

    auto member_yield_col = [&](const std::string &period,
                                const std::string &helicity) -> int {
        auto it = member_yield_col_idx.find(period + "|" + helicity);

        if (it == member_yield_col_idx.end()) return -1;

        return it->second;
    };

    // [FIX] This is the important validity definition. A combined-period member
    // contributes to the combined numerator and denominator only if it has:
    //   acceptance > 0,
    //   acceptance-corrected yield > 0,
    //   finite positive statistical uncertainty,
    //   and positive luminosity for the requested helicity.
    auto member_is_valid_for_combination =
        [&](const std::string &period,
            const std::string &helicity,
            const std::vector<std::string> &fields) -> bool {

            if (!period_has_acceptance(period, fields)) return false;

            const int iy = member_yield_col(period, helicity);

            if (iy < 0 || iy >= (int)fields.size()) return false;

            const Triple Y = parse_tuple3(fields[iy]);

            if (!std::isfinite(Y.value) || Y.value <= 0.0) return false;
            if (!std::isfinite(Y.stat)  || Y.stat  <= 0.0) return false;

            const Triple L = member_period_lumi(period);
            const double lumi_val = lumi_component_for_helicity(L, helicity);

            if (!std::isfinite(lumi_val) || lumi_val <= 0.0) return false;

            return true;
        };

    // [FIX] Combined luminosity now uses the same validity mask as the combined
    // yield numerator. The old code used only acceptance > 0.
    auto lumi_for_label_row = [&](const std::string &L,
                                  const std::string &helicity,
                                  const std::vector<std::string> &fields) -> Triple {
        if (!is_combined_label(L)) {
            auto it = lumi_map.find(L);

            if (it == lumi_map.end()) return Triple{0.0, 0.0, 0.0};

            return it->second;
        }

        Triple out{0.0, 0.0, 0.0};

        for (const auto &m : combined_members_for(L)) {
            if (!member_is_valid_for_combination(m, helicity, fields)) continue;

            auto itLm = lumi_map.find(m);

            if (itLm == lumi_map.end()) continue;

            out.value += itLm->second.value;
            out.stat  += itLm->second.stat;
            out.sys   += itLm->second.sys;
        }

        return out;
    };

    // [FIX] Luminosity columns are scalar/unpolarized columns, so they are written
    // using the unpolarized validity mask.
    auto write_lumi_columns_for_row = [&](std::vector<std::string> &fields) {
        for (const auto &L : labels) {
            const std::string name = lumi_col_for_label(L);

            if (name.empty()) continue;

            auto itc = lumi_col_idx.find(name);

            if (itc == lumi_col_idx.end()) continue;

            Triple lum = lumi_for_label_row(L, "unpol", fields);
            fields[itc->second] = tuple3_to_cell(lum.value, lum.stat, lum.sys);
        }
    };

    // [FIX] Recompute combined acceptance-corrected yields from valid member
    // periods. This intentionally overwrites the upstream combined yield cell
    // for Fa18, Sp18, and 10.6 GeV.
    auto recompute_combined_yield_for_label =
        [&](const std::string &L,
            const std::string &helicity,
            std::vector<std::string> &fields) -> bool {

            if (!is_combined_label(L)) return true;

            auto it_out = colmap.find(L + "|" + helicity);

            if (it_out == colmap.end()) return true;

            const int combined_yield_idx = it_out->second.yield_idx;

            std::vector<Triple> valid_terms;
            valid_terms.reserve(combined_members_for(L).size());

            for (const auto &m : combined_members_for(L)) {
                if (!member_is_valid_for_combination(m, helicity, fields)) continue;

                const int iy = member_yield_col(m, helicity);

                if (iy < 0 || iy >= (int)fields.size()) continue;

                valid_terms.push_back(parse_tuple3(fields[iy]));
            }

            if (valid_terms.empty()) {
                fields[combined_yield_idx] = "";
                return false;
            }

            const Triple combined = add_triples_quadrature_errors(valid_terms);
            fields[combined_yield_idx] = tuple3_to_cell(combined.value,
                                                        combined.stat,
                                                        combined.sys);
            return true;
        };

    std::vector<std::string> out_lines;
    out_lines.reserve(lines.size());
    out_lines.push_back(lines[0]);

    std::cout << "[cross_sections] compute_cross_sections: data rows = "
              << n_data_rows << "\n";

    std::cout << "[cross_sections] NOTE: no imports/efficiency.json correction is applied here. "
              << "Current-efficiency and eppi0 normalization corrections are already upstream.\n";

    std::cout << "[cross_sections] NOTE: accumulated charge is converted to physical integrated "
              << "luminosity with 1.316875 nb^-1/nC before the cross section is formed; "
              << "cross sections are in nb/(GeV^4 deg).\n";

    std::cout << "[cross_sections] NOTE: combined-label luminosities are row-dependent and "
              << "now gated by the same positive-yield validity mask used for the combined-yield numerator.\n";

    std::cout << "[cross_sections] NOTE: combined-label acceptance-corrected yields are recomputed "
              << "inside cross_sections.cpp from valid member periods only.\n";

    std::cout << "[cross_sections] NOTE: Frad/Fbin are imported from Lee's CSV for both energies; "
              << "bin_volume is read from the pass-2 CSV phase-space columns filled by bin_volume.cpp.\n";

    int next_pct = 10;

    for (size_t row = 1; row < lines.size(); ++row) {
        if (lines[row].empty()) {
            out_lines.push_back(lines[row]);
            continue;
        }

        if (n_data_rows > 0) {
            const int pct = (int)std::floor(100.0 * (double)row / (double)n_data_rows);

            if (pct >= next_pct) {
                std::cout << "[cross_sections] compute_cross_sections: ~"
                          << next_pct << "% rows processed\n";
                next_pct += 10;
            }
        }

        std::vector<std::string> fields = split_csv_line(lines[row]);

        if (fields.size() != header.size()) {
            std::cerr << "[cross_sections] WARNING: row " << row
                      << " has " << fields.size() << " fields; expected "
                      << header.size() << ". Copying unchanged.\n";
            out_lines.push_back(lines[row]);
            continue;
        }

        const LeeCorrections *lee_row = nullptr;

        try {
            lee_row = &lee_for_pass2_row(lee_corrections,
                                         fields,
                                         c_pass2_bin_index,
                                         row);
        } catch (const std::exception &e) {
            std::cerr << e.what() << "\n";
            return false;
        }

        const Triple frad = lee_row->frad;
        const Triple fbin = lee_row->fbin;

        // Lee provides the Frad/Fbin values; write those into both energy columns.
        fields[c_frad_106] = tuple3_to_cell(frad.value, frad.stat, frad.sys);
        fields[c_frad_102] = tuple3_to_cell(frad.value, frad.stat, frad.sys);
        fields[c_fbin_106] = tuple3_to_cell(fbin.value, fbin.stat, fbin.sys);
        fields[c_fbin_102] = tuple3_to_cell(fbin.value, fbin.stat, fbin.sys);

        // [FIX] Recompute combined yields before luminosity columns and cross
        // sections are written, so numerator and denominator share the same mask.
        for (const auto &L : labels) {
            if (!is_combined_label(L)) continue;

            for (const auto &h : cross_section_helicities_for_label(L)) {
                (void)recompute_combined_yield_for_label(L, h, fields);
            }
        }

        write_lumi_columns_for_row(fields);

        // The phase-space and diagnostic cubic volumes are identical for all
        // labels/helicities at a given beam energy on this CSV row. Load and
        // validate each energy lazily at most once instead of reparsing the
        // same tuple for every cross-section column.
        struct RowVolumeCache {
            bool loaded = false;
            Triple allowed{0.0, 0.0, 0.0};
        };

        RowVolumeCache volume_106;
        RowVolumeCache volume_102;

        auto row_volume_for_energy = [&](bool use_10p2,
                                         const std::string &context) -> const Triple& {
            RowVolumeCache &cache = use_10p2 ? volume_102 : volume_106;
            if (cache.loaded) return cache.allowed;

            const int c_vbin = use_10p2 ? c_vbin_102 : c_vbin_106;
            const int c_cubic_vbin = use_10p2 ? c_cubic_vbin_102 : c_cubic_vbin_106;
            const std::string energy_tag = use_10p2 ? "10.2 GeV" : "10.6 GeV";

            cache.allowed = parse_required_triple_cell(
                fields,
                c_vbin,
                "bin_volume, " + energy_tag,
                context
            );

            // Preserve the original schema validation. The cubic volume is
            // diagnostic only and is intentionally not used in sigma.
            (void)parse_required_triple_cell(
                fields,
                c_cubic_vbin,
                "cubic bin_volume, " + energy_tag,
                context
            );

            cache.loaded = true;
            return cache.allowed;
        };

        for (const auto &L : labels) {
            const bool use_10p2 = (L == "Sp19 Inb" || L == "10.2 GeV");

            const Triple &Frad = frad;
            const Triple &Fbin = fbin;

            if (Frad.value <= 0.0 || Fbin.value <= 0.0) continue;

            for (const auto &h : cross_section_helicities_for_label(L)) {
                auto it = colmap.find(L + "|" + h);

                if (it == colmap.end()) continue;

                const int iy = it->second.yield_idx;
                const int ix = it->second.xs_idx;

                // Clear the destination first so rerunning this module cannot
                // leave stale positive cross sections in bins that are no longer valid.
                fields[ix] = "";

                const Triple Y = parse_tuple3(fields[iy]);

                if (Y.value <= 0.0 || !std::isfinite(Y.value)) continue;

                const std::string volume_context =
                    "pass-2 CSV row " + std::to_string(row) +
                    ", label " + L + ", helicity " + h;
                const Triple &Vbin = row_volume_for_energy(use_10p2, volume_context);

                const Triple Lumi = lumi_for_label_row(L, h, fields);
                const double lumi_val = lumi_component_for_helicity(Lumi, h);

                if (lumi_val <= 0.0 || !std::isfinite(lumi_val)) continue;

                const double luminosity_nb_inv =
                    integrated_luminosity_nb_inv(lumi_val);
                const double denom = luminosity_nb_inv * Vbin.value;

                if (denom <= 0.0 || !std::isfinite(denom)) continue;

                const double sigma = Y.value * Frad.value * Fbin.value / denom;

                if (!std::isfinite(sigma) || sigma <= 0.0) continue;

                double stat_rel2 = 0.0;
                stat_rel2 += ratio_rel2(Y.value, Y.stat);
                stat_rel2 += ratio_rel2(Frad.value, Frad.stat);
                stat_rel2 += ratio_rel2(Fbin.value, Fbin.stat);
                stat_rel2 += ratio_rel2(Vbin.value, Vbin.stat);

                double sys_rel2 = 0.0;
                sys_rel2 += ratio_rel2(Y.value, Y.sys);
                sys_rel2 += ratio_rel2(Frad.value, Frad.sys);
                sys_rel2 += ratio_rel2(Fbin.value, Fbin.sys);
                sys_rel2 += ratio_rel2(Vbin.value, Vbin.sys);

                const double sigma_stat =
                    (stat_rel2 > 0.0) ? sigma * std::sqrt(stat_rel2) : 0.0;

                const double sigma_sys =
                    (sys_rel2 > 0.0) ? sigma * std::sqrt(sys_rel2) : 0.0;

                fields[ix] = tuple3_to_cell(sigma, sigma_stat, sigma_sys);
            }
        }

        out_lines.push_back(join_csv_line(fields));
    }

    fs::path csv_path(csv_main);
    fs::path tmp_path = csv_path;
    tmp_path += ".tmp";

    std::ofstream ofs(tmp_path);

    if (!ofs) {
        std::cerr << "[cross_sections] ERROR: cannot open "
                  << tmp_path << " for writing.\n";
        return false;
    }

    for (const auto &lout : out_lines) {
        ofs << lout << "\n";
    }

    ofs.close();

    if (!ofs) {
        std::cerr << "[cross_sections] ERROR: failed while writing "
                  << tmp_path << "\n";
        return false;
    }

    std::error_code ec;
    fs::rename(tmp_path, csv_path, ec);

    if (ec) {
        fs::remove(csv_path, ec);
        ec.clear();
        fs::rename(tmp_path, csv_path, ec);

        if (ec) {
            std::cerr << "[cross_sections] ERROR: failed to replace "
                      << csv_main << " with " << tmp_path << ": "
                      << ec.message() << "\n";
            return false;
        }
    }

    std::cout << "[cross_sections] Updated CSV with luminosities and cross sections: "
              << csv_main << "\n";

    return true;
}


// -----------------------------------------------------------------------------
// Analysis-note outputs for the cross-section construction
// -----------------------------------------------------------------------------

static double read_required_scalar_or_tuple_value(
    const std::vector<std::string> &fields,
    int idx) {

    if (idx < 0 || idx >= (int)fields.size()) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    const std::string cell = trim(unquote(fields[idx]));
    if (cell.empty()) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    if (cell.front() == '(') {
        return parse_tuple3(fields[idx]).value;
    }

    return std::atof(cell.c_str());
}

bool write_cross_section_analysis_note_outputs(
    const std::string &csv_main,
    const LumiMap &lumi_map,
    const std::string &out_dir) {

    std::error_code ec;
    fs::create_directories(out_dir, ec);

    if (ec) {
        std::cerr << "[cross_sections] ERROR: cannot create analysis-note directory "
                  << out_dir << ": " << ec.message() << "\n";
        return false;
    }

    {
        std::ofstream o(fs::path(out_dir) / "cross_section_luminosity_summary.csv");
        if (!o) return false;

        o << "run period,beam energy (GeV),selected runs,"
             "accumulated charge (mC),integrated luminosity (pb^-1),"
             "integrated luminosity (fb^-1)\n";

        const std::map<std::string, double> energy = {
            {"Sp18 Inb", 10.594}, {"Sp18 Out", 10.594},
            {"Fa18 Inb", 10.604}, {"Fa18 Out", 10.604},
            {"Sp19 Inb", 10.200}
        };

        const auto &selection = final_pass2_run_selection();

        for (const auto &period : std::vector<std::string>{
                 "Sp18 Inb", "Sp18 Out", "Fa18 Inb", "Fa18 Out", "Sp19 Inb"}) {

            const auto itL = lumi_map.find(period);
            const auto itR = selection.find(period);
            if (itL == lumi_map.end() || itR == selection.end()) continue;

            const double charge_mC = itL->second.value / 1.0e6;
            const double lint_pb = integrated_luminosity_pb_inv(itL->second.value);

            o << period << ","
              << std::fixed << std::setprecision(3) << energy.at(period) << ","
              << itR->second.size() << ","
              << std::setprecision(6) << charge_mC << ","
              << std::setprecision(6) << lint_pb << ","
              << std::setprecision(6) << lint_pb / 1000.0 << "\n";
        }
    }

    {
        std::ofstream o(fs::path(out_dir) / "cross_section_input_summary.csv");
        if (!o) return false;

        o << "symbol,quantity,production source or definition\n";
        o << "N_corr,acceptance/unfolding-corrected DVCS yield,"
             "acceptance corrected yield ep->epg exp <label> unpol\n";
        o << "Q,selected accumulated Faraday-cup charge,"
             "imports/integrated_luminosity/global.csv column 2 summed over final Pass-2 runs\n";
        o << "L_int,integrated luminosity,Q[nC] * 1.316875 nb^-1/nC\n";
        o << "V_bin,physical four-dimensional bin volume,"
             "bin_volume <beam energy> from bin_volume.cpp\n";
        o << "F_rad,radiative correction factor,"
             "Frad imported by bin index from imports/all_bin_v3.csv\n";
        o << "F_bin,bin-centering correction factor,"
             "Fbin imported by bin index from imports/all_bin_v3.csv\n";
        o << "sigma,final four-fold cross section,"
             "N_corr * F_rad * F_bin / (L_int * V_bin)\n";
    }

    std::ifstream ifs(csv_main);
    if (!ifs) return false;

    std::vector<std::string> lines;
    std::string line;
    while (std::getline(ifs, line)) lines.push_back(line);
    if (lines.size() < 2) return false;

    const std::vector<std::string> header = split_csv_line(lines[0]);

    const int c_xbmin = find_col_optional(header, "xBmin");
    const int c_xbmax = find_col_optional(header, "xBmax");
    const int c_q2min = find_col_optional(header, "Q2min");
    const int c_q2max = find_col_optional(header, "Q2max");
    const int c_tmin  = find_col_optional(header, "t_abs_min");
    const int c_tmax  = find_col_optional(header, "t_abs_max");
    const int c_phimin = find_col_optional(header, "phimin");
    const int c_phimax = find_col_optional(header, "phimax");
    const int c_phimean = find_col_optional(header, "phi_mean");

    const std::string combined_yield_column =
        "acceptance corrected yield, ep->epg, exp, "
        + yield_label_for("10.6 GeV")
        + ", unpol";

    const int c_yield =
        find_col_optional(header, combined_yield_column);
    const int c_xs = find_col_optional(
        header, "cross sections, ep->epg, exp, 10.6 GeV, unpol");
    const int c_vbin = find_col_optional(header, "bin_volume, 10.6 GeV");
    const int c_frad = find_col_optional(header, "Frad, 10.6 GeV");
    const int c_fbin = find_col_optional(header, "Fbin, 10.6 GeV");
    const int c_lumi = find_col_optional(
        header, "integrated luminosity, 10.6 GeV (nC)");

    const std::vector<int> required = {
        c_xbmin,c_xbmax,c_q2min,c_q2max,c_tmin,c_tmax,
        c_phimin,c_phimax,c_yield,c_xs,c_vbin,c_frad,c_fbin,c_lumi
    };
    if (c_yield < 0) {
        std::cerr << "[cross_sections] ERROR: missing required analysis-note yield column: "
                  << combined_yield_column << "\n";
        return false;
    }

    for (int idx : required) {
        if (idx < 0) {
            std::cerr << "[cross_sections] ERROR: another required column is missing for "
                         "the cross-section analysis-note example.\n";
            return false;
        }
    }

    struct ExampleRow {
        double xbmin=0, xbmax=0, q2min=0, q2max=0, tmin=0, tmax=0;
        double phimin=0, phimax=0, phi=0;
        Triple yield{0,0,0}, vbin{0,0,0}, frad{0,0,0}, fbin{0,0,0};
        Triple lumi_charge{0,0,0}, xs{0,0,0};
    };

    using KinKey = std::tuple<double,double,double,double,double,double>;
    std::map<KinKey, std::vector<ExampleRow>> groups;

    for (size_t i = 1; i < lines.size(); ++i) {
        if (lines[i].empty()) continue;
        const std::vector<std::string> f = split_csv_line(lines[i]);
        if (f.size() != header.size()) continue;

        const Triple Y = parse_tuple3(f[c_yield]);
        const Triple V = parse_tuple3(f[c_vbin]);
        const Triple R = parse_tuple3(f[c_frad]);
        const Triple B = parse_tuple3(f[c_fbin]);
        const Triple L = parse_tuple3(f[c_lumi]);
        const Triple X = parse_tuple3(f[c_xs]);

        if (!(Y.value > 0.0) || !(V.value > 0.0) ||
            !(R.value > 0.0) || !(B.value > 0.0) ||
            !(L.value > 0.0) || !(X.value > 0.0)) continue;

        ExampleRow r;
        r.xbmin = read_required_scalar_or_tuple_value(f,c_xbmin);
        r.xbmax = read_required_scalar_or_tuple_value(f,c_xbmax);
        r.q2min = read_required_scalar_or_tuple_value(f,c_q2min);
        r.q2max = read_required_scalar_or_tuple_value(f,c_q2max);
        r.tmin = read_required_scalar_or_tuple_value(f,c_tmin);
        r.tmax = read_required_scalar_or_tuple_value(f,c_tmax);
        r.phimin = read_required_scalar_or_tuple_value(f,c_phimin);
        r.phimax = read_required_scalar_or_tuple_value(f,c_phimax);
        r.phi = (c_phimean >= 0)
            ? read_required_scalar_or_tuple_value(f,c_phimean)
            : 0.5*(r.phimin+r.phimax);
        r.yield=Y; r.vbin=V; r.frad=R; r.fbin=B; r.lumi_charge=L; r.xs=X;

        groups[KinKey{r.xbmin,r.xbmax,r.q2min,r.q2max,r.tmin,r.tmax}].push_back(r);
    }

    if (groups.empty()) {
        std::cerr << "[cross_sections] WARNING: no complete 10.6-GeV group found "
                     "for analysis-note example.\n";
        return true;
    }

    // ---------------------------------------------------------------------
    // Build period-level raw-count inputs needed for the first stage of the
    // illustrative correction chain.  The combined 10.6-GeV raw count in a
    // row uses the same member-period validity mask as the combined
    // acceptance-corrected yield.
    // ---------------------------------------------------------------------
    const std::vector<std::string> periods_10p6 = {
        "Fa18 Inb", "Fa18 Out", "Sp18 Inb", "Sp18 Out"
    };
    const std::vector<std::string> topologies = {
        "(FD, FD)", "(CD, FD)", "(CD, FT)"
    };

    std::map<std::string,int> c_period_acceptance;
    std::map<std::string,std::vector<int>> c_period_raw;

    for(const auto &period:periods_10p6){
        c_period_acceptance[period] =
            find_col_optional(
                header,
                "acceptance corrected yield, ep->epg, exp, "
                + period + ", unpol"
            );

        std::vector<int> cols;
        for(const auto &topology:topologies){
            cols.push_back(
                find_col_optional(
                    header,
                    "raw yield, ep->epg, "
                    + topology + ", exp, "
                    + period + ", unpol"
                )
            );
        }
        c_period_raw[period]=cols;
    }

    for(const auto &period:periods_10p6){
        if(c_period_acceptance[period]<0){
            std::cerr
                << "[cross_sections] ERROR: missing period acceptance column for "
                << period << " in analysis-note correction-chain output.\n";
            return false;
        }
        for(int idx:c_period_raw[period]){
            if(idx<0){
                std::cerr
                    << "[cross_sections] ERROR: missing raw-yield topology column for "
                    << period << " in analysis-note correction-chain output.\n";
                return false;
            }
        }
    }

    struct GroupCandidate {
        KinKey key;
        std::vector<ExampleRow> rows;
        double xb=0.0, q2=0.0, tt=0.0;
    };

    std::vector<GroupCandidate> candidates;
    candidates.reserve(groups.size());

    double xb_lo=std::numeric_limits<double>::infinity();
    double xb_hi=-std::numeric_limits<double>::infinity();
    double q2_lo=std::numeric_limits<double>::infinity();
    double q2_hi=-std::numeric_limits<double>::infinity();
    double t_lo=std::numeric_limits<double>::infinity();
    double t_hi=-std::numeric_limits<double>::infinity();

    std::size_t max_phi_points=0;

    for(const auto &kv:groups){
        GroupCandidate cnd;
        cnd.key=kv.first;
        cnd.rows=kv.second;
        cnd.xb=0.5*(std::get<0>(kv.first)+std::get<1>(kv.first));
        cnd.q2=0.5*(std::get<2>(kv.first)+std::get<3>(kv.first));
        cnd.tt=0.5*(std::get<4>(kv.first)+std::get<5>(kv.first));

        std::sort(
            cnd.rows.begin(),cnd.rows.end(),
            [](const ExampleRow&a,const ExampleRow&b){return a.phi<b.phi;}
        );

        max_phi_points=std::max(max_phi_points,cnd.rows.size());
        xb_lo=std::min(xb_lo,cnd.xb);
        xb_hi=std::max(xb_hi,cnd.xb);
        q2_lo=std::min(q2_lo,cnd.q2);
        q2_hi=std::max(q2_hi,cnd.q2);
        t_lo=std::min(t_lo,cnd.tt);
        t_hi=std::max(t_hi,cnd.tt);

        candidates.push_back(std::move(cnd));
    }

    auto unit=[](double x,double lo,double hi){
        return (hi>lo)?(x-lo)/(hi-lo):0.5;
    };

    // Four deliberately separated locations in the occupied 10.6-GeV phase
    // space.  The nearest available cell is used, with strong preference for
    // complete/broad phi coverage.
    const std::array<std::array<double,3>,4> targets = {{
        {{0.12,0.12,0.15}},
        {{0.36,0.34,0.32}},
        {{0.62,0.62,0.58}},
        {{0.88,0.88,0.82}}
    }};

    std::vector<GroupCandidate> examples;
    std::set<KinKey> used_keys;

    for(const auto &target:targets){
        bool have=false;
        double best_score=std::numeric_limits<double>::infinity();
        GroupCandidate best_candidate;

        for(const auto &cnd:candidates){
            if(used_keys.count(cnd.key)) continue;

            const double ux=unit(cnd.xb,xb_lo,xb_hi);
            const double uq=unit(cnd.q2,q2_lo,q2_hi);
            const double ut=unit(cnd.tt,t_lo,t_hi);

            const double distance =
                std::pow(ux-target[0],2)
                +std::pow(uq-target[1],2)
                +std::pow(ut-target[2],2);

            const double coverage_penalty =
                (max_phi_points>0)
                ? 2.0*(1.0-static_cast<double>(cnd.rows.size())
                            /static_cast<double>(max_phi_points))
                : 0.0;

            const double score=distance+coverage_penalty;

            if(score<best_score){
                best_score=score;
                best_candidate=cnd;
                have=true;
            }
        }

        if(have){
            used_keys.insert(best_candidate.key);
            examples.push_back(std::move(best_candidate));
        }
    }

    if(examples.size()<4){
        std::sort(
            candidates.begin(),candidates.end(),
            [](const GroupCandidate&a,const GroupCandidate&b){
                return a.rows.size()>b.rows.size();
            }
        );

        for(const auto &cnd:candidates){
            if(examples.size()>=4) break;
            if(used_keys.count(cnd.key)) continue;
            used_keys.insert(cnd.key);
            examples.push_back(cnd);
        }
    }

    // ---------------------------------------------------------------------
    // Re-read the CSV into a row lookup so the raw selected counts can be
    // reconstructed for the same rows used by the representative examples.
    // ---------------------------------------------------------------------
    struct RawStageRow {
        double phi=0.0;
        double raw_total=0.0;
        double raw_stat=0.0;
    };

    std::map<std::tuple<double,double,double,double,double,double,double,double>,
             RawStageRow> raw_stage_lookup;

    for(size_t i=1;i<lines.size();++i){
        if(lines[i].empty()) continue;

        const std::vector<std::string> f=split_csv_line(lines[i]);
        if(f.size()!=header.size()) continue;

        const double xbmin=read_required_scalar_or_tuple_value(f,c_xbmin);
        const double xbmax=read_required_scalar_or_tuple_value(f,c_xbmax);
        const double q2min=read_required_scalar_or_tuple_value(f,c_q2min);
        const double q2max=read_required_scalar_or_tuple_value(f,c_q2max);
        const double ttmin=read_required_scalar_or_tuple_value(f,c_tmin);
        const double ttmax=read_required_scalar_or_tuple_value(f,c_tmax);
        const double phimin=read_required_scalar_or_tuple_value(f,c_phimin);
        const double phimax=read_required_scalar_or_tuple_value(f,c_phimax);

        double raw_total=0.0;
        double raw_var=0.0;

        for(const auto &period:periods_10p6){
            const Triple period_acceptance =
                parse_tuple3(f[c_period_acceptance[period]]);

            // Match the combined-yield validity rule: only member periods that
            // contributed a valid acceptance-corrected yield contribute their
            // raw selected counts here.
            if(!(period_acceptance.value>0.0)) continue;

            for(int col:c_period_raw[period]){
                const Triple raw_count = parse_tuple3(f[col]);

                if(std::isfinite(raw_count.value) && raw_count.value>=0.0){
                    raw_total += raw_count.value;

                    if(std::isfinite(raw_count.stat) && raw_count.stat>=0.0){
                        raw_var += raw_count.stat * raw_count.stat;
                    }
                }
            }
        }

        RawStageRow rr;
        rr.phi=(c_phimean>=0)
            ? read_required_scalar_or_tuple_value(f,c_phimean)
            : 0.5*(phimin+phimax);
        rr.raw_total=raw_total;
        rr.raw_stat=std::sqrt(std::max(0.0,raw_var));

        raw_stage_lookup[
            std::make_tuple(
                xbmin,xbmax,q2min,q2max,ttmin,ttmax,phimin,phimax
            )
        ]=rr;
    }

    // ---------------------------------------------------------------------
    // Numerical table for all five displayed stages.
    // Every stage is expressed in cross-section units using the same
    // L_int*V_bin denominator.  This makes the log-scale overlay meaningful.
    // ---------------------------------------------------------------------
    {
        std::ofstream o(
            fs::path(out_dir)/"cross_section_assembly_example.csv"
        );
        if(!o) return false;

        o<<"example,xBmin,xBmax,Q2min,Q2max,t_abs_min,t_abs_max,"
           "phi_min,phi_max,phi_mean,"
           "total_selected_counts,total_selected_counts_stat,"
           "acceptance_corrected_yield,acceptance_corrected_yield_stat,"
           "charge_10p6_nC,Lint_10p6_nb^-1,V_bin,F_rad,F_bin,"
           "stage_total_counts_nb_per_GeV4_deg,"
           "stage_acceptance_corrected_nb_per_GeV4_deg,"
           "stage_after_Frad_nb_per_GeV4_deg,"
           "stage_after_Fbin_nb_per_GeV4_deg,"
           "final_stored_sigma_nb_per_GeV4_deg,final_stored_sigma_stat\n";

        o<<std::setprecision(12);

        for(std::size_t iex=0;iex<examples.size();++iex){
            const char example_label=static_cast<char>('A'+iex);

            for(const auto&r:examples[iex].rows){
                const auto raw_key=
                    std::make_tuple(
                        r.xbmin,r.xbmax,r.q2min,r.q2max,r.tmin,r.tmax,
                        r.phimin,r.phimax
                    );

                const auto it_raw=raw_stage_lookup.find(raw_key);
                if(it_raw==raw_stage_lookup.end()) continue;

                const double lint=
                    integrated_luminosity_nb_inv(r.lumi_charge.value);
                const double denom=lint*r.vbin.value;

                if(!(denom>0.0)) continue;

                const double stage_raw=
                    it_raw->second.raw_total/denom;
                const double stage_acceptance=
                    r.yield.value/denom;
                const double stage_rad=
                    stage_acceptance*r.frad.value;
                const double stage_bin=
                    stage_rad*r.fbin.value;

                o<<example_label<<","
                 <<r.xbmin<<","<<r.xbmax<<","
                 <<r.q2min<<","<<r.q2max<<","
                 <<r.tmin<<","<<r.tmax<<","
                 <<r.phimin<<","<<r.phimax<<","<<r.phi<<","
                 <<it_raw->second.raw_total<<","
                 <<it_raw->second.raw_stat<<","
                 <<r.yield.value<<","<<r.yield.stat<<","
                 <<r.lumi_charge.value<<","<<lint<<","
                 <<r.vbin.value<<","<<r.frad.value<<","<<r.fbin.value<<","
                 <<stage_raw<<","
                 <<stage_acceptance<<","
                 <<stage_rad<<","
                 <<stage_bin<<","
                 <<r.xs.value<<","<<r.xs.stat<<"\n";
            }
        }
    }

    // ---------------------------------------------------------------------
    // 2x2 log-scale illustration of the full extraction chain.
    // ---------------------------------------------------------------------
    {
        TCanvas c(
            "c_cross_section_chain_note","",
            1500,1100
        );
        c.Divide(2,2,0.002,0.002);

        for(std::size_t iex=0;iex<examples.size() && iex<4;++iex){
            c.cd(static_cast<int>(iex)+1);

            gPad->SetLeftMargin((iex%2==0)?0.145:0.115);
            gPad->SetRightMargin(0.035);
            gPad->SetBottomMargin((iex>=2)?0.145:0.115);
            gPad->SetTopMargin(0.14);
            gPad->SetTicks(1,1);
            gPad->SetLogy();

            TGraphErrors g_raw;
            TGraphErrors g_acc;
            TGraphErrors g_rad;
            TGraphErrors g_bin;
            TGraphErrors g_final;

            double ymin=std::numeric_limits<double>::infinity();
            double ymax=0.0;
            int ip=0;

            for(const auto&r:examples[iex].rows){
                const auto raw_key=
                    std::make_tuple(
                        r.xbmin,r.xbmax,r.q2min,r.q2max,r.tmin,r.tmax,
                        r.phimin,r.phimax
                    );

                const auto it_raw=raw_stage_lookup.find(raw_key);
                if(it_raw==raw_stage_lookup.end()) continue;

                const double lint=
                    integrated_luminosity_nb_inv(r.lumi_charge.value);
                const double denom=lint*r.vbin.value;
                if(!(denom>0.0)) continue;

                const double raw=
                    it_raw->second.raw_total/denom;
                const double acc=
                    r.yield.value/denom;
                const double after_rad=
                    acc*r.frad.value;
                const double after_bin=
                    after_rad*r.fbin.value;
                const double final=
                    r.xs.value;

                const double raw_err=
                    it_raw->second.raw_stat/denom;

                const double rel_y=
                    (r.yield.value>0.0)
                    ? r.yield.stat/r.yield.value : 0.0;
                const double rel_v=
                    (r.vbin.value>0.0)
                    ? r.vbin.stat/r.vbin.value : 0.0;

                const double acc_err=
                    acc*std::sqrt(
                        rel_y*rel_y+rel_v*rel_v
                    );

                const double rel_r=
                    (r.frad.value>0.0)
                    ? r.frad.stat/r.frad.value : 0.0;
                const double rad_err=
                    after_rad*std::sqrt(
                        rel_y*rel_y+rel_v*rel_v+rel_r*rel_r
                    );

                const double rel_b=
                    (r.fbin.value>0.0)
                    ? r.fbin.stat/r.fbin.value : 0.0;
                const double bin_err=
                    after_bin*std::sqrt(
                        rel_y*rel_y+rel_v*rel_v
                        +rel_r*rel_r+rel_b*rel_b
                    );

                g_raw.SetPoint(ip,r.phi,raw);
                g_raw.SetPointError(ip,0.0,raw_err);

                g_acc.SetPoint(ip,r.phi,acc);
                g_acc.SetPointError(ip,0.0,acc_err);

                g_rad.SetPoint(ip,r.phi,after_rad);
                g_rad.SetPointError(ip,0.0,rad_err);

                g_bin.SetPoint(ip,r.phi,after_bin);
                g_bin.SetPointError(ip,0.0,bin_err);

                g_final.SetPoint(ip,r.phi,final);
                g_final.SetPointError(ip,0.0,r.xs.stat);

                for(double v:{raw,acc,after_rad,after_bin,final}){
                    if(v>0.0 && std::isfinite(v)){
                        ymin=std::min(ymin,v);
                        ymax=std::max(ymax,v);
                    }
                }

                ++ip;
            }

            if(!(ymin>0.0) || !(ymax>ymin)){
                ymin=1.0e-4;
                ymax=1.0;
            }

            const double log_lo=
                std::floor(std::log10(ymin))-0.15;
            const double log_hi=
                std::ceil(std::log10(ymax))+0.15;

            const std::string frame_name=
                "h_cross_section_chain_note_"
                +std::to_string(iex);

            TH1F frame(
                frame_name.c_str(),"",
                100,0.0,360.0
            );
            frame.SetMinimum(std::pow(10.0,log_lo));
            frame.SetMaximum(std::pow(10.0,log_hi));
            frame.GetXaxis()->SetTitle("#phi (deg)");
            frame.GetYaxis()->SetTitle(
                "cross-section-equivalent value  [nb/(GeV^{4} deg)]"
            );
            frame.GetXaxis()->SetTitleSize(0.052);
            frame.GetYaxis()->SetTitleSize(0.045);
            frame.GetXaxis()->SetLabelSize(0.043);
            frame.GetYaxis()->SetLabelSize(0.043);
            frame.GetYaxis()->SetTitleOffset(
                (iex%2==0)?1.50:1.18
            );
            frame.DrawCopy();

            g_raw.SetMarkerStyle(24);
            g_raw.SetMarkerSize(0.88);
            g_raw.SetLineWidth(2);
            g_raw.SetMarkerColor(kGray+2);
            g_raw.SetLineColor(kGray+2);

            g_acc.SetMarkerStyle(25);
            g_acc.SetMarkerSize(0.88);
            g_acc.SetLineWidth(2);
            g_acc.SetMarkerColor(kBlack);
            g_acc.SetLineColor(kBlack);

            g_rad.SetMarkerStyle(26);
            g_rad.SetMarkerSize(0.88);
            g_rad.SetLineWidth(2);
            g_rad.SetMarkerColor(kBlue+1);
            g_rad.SetLineColor(kBlue+1);

            g_bin.SetMarkerStyle(32);
            g_bin.SetMarkerSize(0.92);
            g_bin.SetLineWidth(2);
            g_bin.SetMarkerColor(kMagenta+1);
            g_bin.SetLineColor(kMagenta+1);

            g_final.SetMarkerStyle(20);
            g_final.SetMarkerSize(0.78);
            g_final.SetLineWidth(2);
            g_final.SetMarkerColor(kRed+1);
            g_final.SetLineColor(kRed+1);

            g_raw.DrawClone("PE SAME");
            g_acc.DrawClone("PE SAME");
            g_rad.DrawClone("PE SAME");
            g_bin.DrawClone("PE SAME");
            g_final.DrawClone("PE SAME");

            const KinKey&k=examples[iex].key;
            std::ostringstream kin;
            kin<<std::fixed<<std::setprecision(3)
               <<std::get<0>(k)<<" < x_{B} < "<<std::get<1>(k)
               <<",  "<<std::get<2>(k)<<" < Q^{2} < "
               <<std::get<3>(k)<<" GeV^{2}"
               <<",  "<<std::get<4>(k)<<" < |t| < "
               <<std::get<5>(k)<<" GeV^{2}";

            TLatex latex;
            latex.SetNDC();
            latex.SetTextFont(42);
            latex.SetTextSize(0.036);

            const std::string panel_label=
                std::string("(")
                +static_cast<char>('a'+iex)+")";

            latex.DrawLatex(
                0.15,0.925,
                panel_label.c_str()
            );

            latex.SetTextSize(0.027);
            latex.DrawLatex(
                0.15,0.875,
                kin.str().c_str()
            );

            if(iex==0){
                TLegend leg(
                    0.48,0.54,0.94,0.82
                );
                leg.SetBorderSize(0);
                leg.SetFillStyle(0);
                leg.SetTextSize(0.028);
                leg.AddEntry(
                    &g_raw,
                    "total selected counts / (L_{int} V_{bin})",
                    "pe"
                );
                leg.AddEntry(
                    &g_acc,
                    "acceptance corrected / (L_{int} V_{bin})",
                    "pe"
                );
                leg.AddEntry(
                    &g_rad,
                    "after F_{rad}",
                    "pe"
                );
                leg.AddEntry(
                    &g_bin,
                    "after F_{bin}",
                    "pe"
                );
                leg.AddEntry(
                    &g_final,
                    "final stored cross section",
                    "pe"
                );
                leg.DrawClone();
            }
        }

        c.cd(0);

        TLatex title;
        title.SetNDC();
        title.SetTextFont(42);
        title.SetTextAlign(22);
        title.SetTextSize(0.027);
        title.DrawLatex(
            0.50,0.988,
            "Representative construction of the 10.6 GeV unpolarized DVCS cross section"
        );

        c.SaveAs(
            (
                fs::path(out_dir)
                /"cross_section_correction_chain_example.png"
            ).string().c_str()
        );
    }

    std::cout<<"[cross_sections] Wrote analysis-note outputs to "<<out_dir<<"\n";
    return true;
}

// -----------------------------------------------------------------------------
// Plotting structures and helpers
// -----------------------------------------------------------------------------

struct Point {
    double phi;
    double xs;
    double xs_err;
};

struct BinData {
    std::vector<Point> unpol;
    std::vector<Point> pos;
    std::vector<Point> neg;
    size_t theory_row = 0;
    bool have_theory_row = false;
};

using QTKey = std::pair<Range, Range>;

struct XSGroupByXB {
    std::map<QTKey, BinData> bins;
    int xb_index = -1;
};

struct TheoryCurves {
    std::vector<double> phi_deg;
    std::vector<double> bh_unpol;
    std::vector<double> bh_pos;
    std::vector<double> bh_neg;
    std::vector<double> km_unpol;
    std::vector<double> km_pos;
    std::vector<double> km_neg;
    std::vector<double> vgg_unpol;
    std::vector<double> vgg_pos;
    std::vector<double> vgg_neg;
};

static std::string theory_energy_label_for(const std::string &label) {
    if (label == "Sp19 Inb" || label == "10.2 GeV") return "10.2 GeV";
    return "10.6 GeV";
}

static std::map<size_t, TheoryCurves>
load_theory_for_label(const std::string &label,
                      const std::string &theory_root) {
    std::map<size_t, TheoryCurves> out;

    std::string energy_label = theory_energy_label_for(label);
    fs::path dir  = fs::path(theory_root) / canonical_period_dir(energy_label);
    fs::path file = dir / "xs_phi_all.json";

    if (!fs::exists(file)) {
        std::cerr << "[cross_sections] WARNING: no theory JSON for label \""
                  << label << "\" at " << file.string() << "\n";
        return out;
    }

    std::ifstream ifs(file);

    if (!ifs) {
        std::cerr << "[cross_sections] WARNING: cannot open theory JSON for label \""
                  << label << "\" at " << file.string() << "\n";
        return out;
    }

    json j;

    try {
        ifs >> j;
    } catch (...) {
        std::cerr << "[cross_sections] WARNING: malformed theory JSON for label \""
                  << label << "\" at " << file.string() << "\n";
        return out;
    }

    std::vector<double> phi_deg = j.value("phi_deg", std::vector<double>{});

    if (phi_deg.empty()) {
        std::cerr << "[cross_sections] WARNING: theory JSON for label \""
                  << label << "\" has empty phi_deg.\n";
        return out;
    }

    if (!j.contains("rows") || !j["rows"].is_object()) {
        std::cerr << "[cross_sections] WARNING: theory JSON for label \""
                  << label << "\" has no rows object.\n";
        return out;
    }

    for (auto it = j["rows"].begin(); it != j["rows"].end(); ++it) {
        const std::string row_key = it.key();
        const json &cell = it.value();

        size_t row_index = 0;

        try {
            row_index = (size_t)std::stoul(row_key);
        } catch (...) {
            continue;
        }

        TheoryCurves tc;
        tc.phi_deg   = phi_deg;
        tc.bh_unpol  = cell["BH"].value("unpol", std::vector<double>{});
        tc.bh_pos    = cell["BH"].value("pos",   std::vector<double>{});
        tc.bh_neg    = cell["BH"].value("neg",   std::vector<double>{});
        tc.km_unpol  = cell["KM"].value("unpol", std::vector<double>{});
        tc.km_pos    = cell["KM"].value("pos",   std::vector<double>{});
        tc.km_neg    = cell["KM"].value("neg",   std::vector<double>{});
        tc.vgg_unpol = cell["VGG"].value("unpol", std::vector<double>{});
        tc.vgg_pos   = cell["VGG"].value("pos",   std::vector<double>{});
        tc.vgg_neg   = cell["VGG"].value("neg",   std::vector<double>{});

        if (!tc.phi_deg.empty()) {
            out[row_index] = std::move(tc);
        }
    }

    std::cout << "[cross_sections] Loaded theory for label \"" << label
              << "\" (energy " << energy_label << ") rows=" << out.size()
              << " from " << file.string() << "\n";

    return out;
}

enum class XSecPanelMode {
    All,
    UnpolOnly,
    PosOnly,
    NegOnly
};

static std::pair<double, double> compute_yrange_for_bin(
    const BinData *bin,
    const std::map<size_t, TheoryCurves> &theory,
    XSecPanelMode mode) {

    double ymin = std::numeric_limits<double>::max();
    double ymax = 0.0;

    auto update_from_points = [&](const std::vector<Point> &v) {
        for (const auto &p : v) {
            if (p.xs > 0.0) {
                double ylow  = std::max(1e-12, p.xs - p.xs_err);
                double yhigh = p.xs + p.xs_err;

                if (ylow  > 0.0 && ylow  < ymin) ymin = ylow;
                if (yhigh > ymax) ymax = yhigh;
            }
        }
    };

    auto update_from_curve = [&](const std::vector<double> &ys) {
        for (double y : ys) {
            if (y <= 0.0) continue;

            if (y < ymin) ymin = y;
            if (y > ymax) ymax = y;
        }
    };

    if (bin) {
        if (mode == XSecPanelMode::All || mode == XSecPanelMode::UnpolOnly) {
            update_from_points(bin->unpol);
        }

        if (mode == XSecPanelMode::All || mode == XSecPanelMode::PosOnly) {
            update_from_points(bin->pos);
        }

        if (mode == XSecPanelMode::All || mode == XSecPanelMode::NegOnly) {
            update_from_points(bin->neg);
        }

        if (bin->have_theory_row) {
            auto it_th = theory.find(bin->theory_row);

            if (it_th != theory.end()) {
                const TheoryCurves &tc = it_th->second;

                if (mode == XSecPanelMode::All) {
                    update_from_curve(tc.bh_unpol);
                    update_from_curve(tc.km_unpol);
                    update_from_curve(tc.km_pos);
                    update_from_curve(tc.km_neg);
                    update_from_curve(tc.vgg_unpol);
                    update_from_curve(tc.vgg_pos);
                    update_from_curve(tc.vgg_neg);
                } else if (mode == XSecPanelMode::UnpolOnly) {
                    update_from_curve(tc.bh_unpol);
                    update_from_curve(tc.km_unpol);
                    update_from_curve(tc.vgg_unpol);
                } else if (mode == XSecPanelMode::PosOnly) {
                    update_from_curve(tc.bh_unpol);
                    update_from_curve(tc.km_pos);
                    update_from_curve(tc.vgg_pos);
                } else if (mode == XSecPanelMode::NegOnly) {
                    update_from_curve(tc.bh_unpol);
                    update_from_curve(tc.km_neg);
                    update_from_curve(tc.vgg_neg);
                }
            }
        }
    }

    if (ymax <= 0.0 || !std::isfinite(ymax)) {
        ymin = 1e-4;
        ymax = 1.0;
    } else {
        if (ymin <= 0.0 || !std::isfinite(ymin)) {
            ymin = ymax * 1e-3;
        }

        double logmin = std::pow(10.0, std::floor(std::log10(ymin)));
        double logmax = std::pow(10.0, std::ceil(std::log10(ymax)));

        ymin = std::max(1e-4, logmin);
        ymax = logmax;
    }

    return std::make_pair(ymin, 1.2 * ymax);
}

static TGraphErrors *make_xsec_graph(const std::vector<Point> &v,
                                     int mstyle,
                                     int mcolor) {
    if (v.empty()) return nullptr;

    int N = (int)v.size();

    std::vector<double> x(N);
    std::vector<double> y(N);
    std::vector<double> ex(N);
    std::vector<double> ey(N);

    for (int i = 0; i < N; ++i) {
        x[i]  = v[i].phi;
        y[i]  = v[i].xs;
        ex[i] = 0.0;
        ey[i] = v[i].xs_err;
    }

    TGraphErrors *g = new TGraphErrors(N, x.data(), y.data(), ex.data(), ey.data());
    g->SetMarkerStyle(mstyle);
    g->SetMarkerSize(1.0);
    g->SetLineWidth(2);
    g->SetLineColor(mcolor);
    g->SetMarkerColor(mcolor);

    return g;
}

// -----------------------------------------------------------------------------
// Canvas builder (one xB bin, one view mode)
// -----------------------------------------------------------------------------

static void make_xsec_canvas_for_mode(
    const std::string &label,
    const Range &xb_range,
    const XSGroupByXB &group,
    const std::vector<Range> &q2_slice,
    const std::vector<Range> &t_slice,
    const std::map<size_t, TheoryCurves> &theory,
    const fs::path &outdir,
    int xb_idx_for_name,
    XSecPanelMode mode,
    int ncols,
    int nrows,
    int nPads) {

    const auto &bins_for_xB = group.bins;

    if (bins_for_xB.empty()) return;

    int W = 280 * ncols + 160;
    int H = 260 * nrows + 260;

    if (W < 1200) W = 1200;
    if (H < 900)  H = 900;

    double titleSize       = 0.18;
    double legendTextSize  = 0.11;
    double cellLabelSize   = 0.070;

    if (nPads <= 4) {
        titleSize      = 0.14;
        legendTextSize = 0.09;
        cellLabelSize  = 0.060;
    }

    if (nPads == 1) {
        titleSize      = 0.12;
        legendTextSize = 0.085;
        cellLabelSize  = 0.055;
    }

    titleSize      *= 0.5;
    legendTextSize *= 0.5;

    std::ostringstream cname;
    cname << "c_xsec_";

    if (mode == XSecPanelMode::UnpolOnly) {
        cname << "unpol_";
    } else if (mode == XSecPanelMode::PosOnly) {
        cname << "pos_";
    } else if (mode == XSecPanelMode::NegOnly) {
        cname << "neg_";
    }

    cname << canonical_period_dir(label) << "_xB" << xb_idx_for_name;

    TCanvas *c = new TCanvas(cname.str().c_str(), cname.str().c_str(), W, H);

    TPad *pTop = new TPad("pTop", "pTop", 0.0, 0.78, 1.0, 1.0);
    pTop->SetFillStyle(0);
    pTop->SetBorderSize(0);
    pTop->Draw();

    TPad *pGrid = new TPad("pGrid", "pGrid", 0.0, 0.00, 1.0, 0.78);
    pGrid->SetFillStyle(0);
    pGrid->SetBorderSize(0);
    pGrid->Draw();
    pGrid->cd();
    pGrid->Divide(ncols, nrows, 0.0001, 0.0001);

    pTop->cd();

    TLatex head;
    head.SetNDC();
    head.SetTextAlign(22);
    head.SetTextFont(42);
    head.SetTextSize(titleSize);

    std::ostringstream tit;
    tit << "Cross sections, ep #rightarrow ep#gamma   " << label
        << "   x_{B} in ("
        << std::fixed << std::setprecision(3)
        << xb_range.first << ", " << xb_range.second << ")";

    if (mode == XSecPanelMode::UnpolOnly) {
        tit << "   (unpolarized only)";
    } else if (mode == XSecPanelMode::PosOnly) {
        tit << "   (+ helicity only)";
    } else if (mode == XSecPanelMode::NegOnly) {
        tit << "   (- helicity only)";
    }

    head.DrawLatex(0.5, 0.86, tit.str().c_str());

    TGraphErrors dummy_unpol;
    TGraphErrors dummy_pos;
    TGraphErrors dummy_neg;

    dummy_unpol.SetMarkerStyle(20);
    dummy_unpol.SetMarkerSize(1.0);
    dummy_unpol.SetLineWidth(2);
    dummy_unpol.SetMarkerColor(kBlack);
    dummy_unpol.SetLineColor(kBlack);

    dummy_pos.SetMarkerStyle(24);
    dummy_pos.SetMarkerSize(1.0);
    dummy_pos.SetLineWidth(2);
    dummy_pos.SetMarkerColor(kRed + 1);
    dummy_pos.SetLineColor(kRed + 1);

    dummy_neg.SetMarkerStyle(25);
    dummy_neg.SetMarkerSize(1.0);
    dummy_neg.SetLineWidth(2);
    dummy_neg.SetMarkerColor(kBlue + 1);
    dummy_neg.SetLineColor(kBlue + 1);

    TGraph dummy_bh;
    TGraph dummy_km_unpol;
    TGraph dummy_km_pos;
    TGraph dummy_km_neg;
    TGraph dummy_vgg_unpol;
    TGraph dummy_vgg_pos;
    TGraph dummy_vgg_neg;

    dummy_bh.SetLineWidth(2);
    dummy_bh.SetLineStyle(2);
    dummy_bh.SetLineColor(kGreen + 2);

    dummy_km_unpol.SetLineWidth(2);
    dummy_km_unpol.SetLineStyle(1);
    dummy_km_unpol.SetLineColor(kMagenta + 1);

    dummy_km_pos.SetLineWidth(2);
    dummy_km_pos.SetLineStyle(2);
    dummy_km_pos.SetLineColor(kMagenta + 1);

    dummy_km_neg.SetLineWidth(2);
    dummy_km_neg.SetLineStyle(3);
    dummy_km_neg.SetLineColor(kMagenta + 1);

    dummy_vgg_unpol.SetLineWidth(2);
    dummy_vgg_unpol.SetLineStyle(1);
    dummy_vgg_unpol.SetLineColor(kOrange + 7);

    dummy_vgg_pos.SetLineWidth(2);
    dummy_vgg_pos.SetLineStyle(2);
    dummy_vgg_pos.SetLineColor(kOrange + 7);

    dummy_vgg_neg.SetLineWidth(2);
    dummy_vgg_neg.SetLineStyle(3);
    dummy_vgg_neg.SetLineColor(kOrange + 7);

    TLegend *legData = new TLegend(0.02, 0.05, 0.32, 0.80);
    legData->SetBorderSize(1);
    legData->SetLineColor(kBlack);
    legData->SetFillColor(kWhite);
    legData->SetFillStyle(1001);
    legData->SetTextFont(42);
    legData->SetTextSize(legendTextSize);

    TLegend *legKM = new TLegend(0.35, 0.05, 0.65, 0.80);
    legKM->SetBorderSize(1);
    legKM->SetLineColor(kBlack);
    legKM->SetFillColor(kWhite);
    legKM->SetFillStyle(1001);
    legKM->SetTextFont(42);
    legKM->SetTextSize(legendTextSize);

    TLegend *legVGG = new TLegend(0.68, 0.05, 0.98, 0.80);
    legVGG->SetBorderSize(1);
    legVGG->SetLineColor(kBlack);
    legVGG->SetFillColor(kWhite);
    legVGG->SetFillStyle(1001);
    legVGG->SetTextFont(42);
    legVGG->SetTextSize(legendTextSize);

    if (mode == XSecPanelMode::All) {
        legData->AddEntry(&dummy_unpol, "data unpolarized", "lep");
        legData->AddEntry(&dummy_pos,   "data + helicity", "lep");
        legData->AddEntry(&dummy_neg,   "data - helicity", "lep");

        legKM->AddEntry(&dummy_bh,       "BH unpolarized", "l");
        legKM->AddEntry(&dummy_km_unpol, "KM unpolarized", "l");
        legKM->AddEntry(&dummy_km_pos,   "KM + helicity",  "l");
        legKM->AddEntry(&dummy_km_neg,   "KM - helicity",  "l");

        legVGG->AddEntry(&dummy_vgg_unpol, "VGG unpolarized", "l");
        legVGG->AddEntry(&dummy_vgg_pos,   "VGG + helicity",  "l");
        legVGG->AddEntry(&dummy_vgg_neg,   "VGG - helicity",  "l");
    } else if (mode == XSecPanelMode::UnpolOnly) {
        legData->AddEntry(&dummy_unpol, "data unpolarized", "lep");
        legKM->AddEntry(&dummy_bh,       "BH unpolarized", "l");
        legKM->AddEntry(&dummy_km_unpol, "KM unpolarized", "l");
        legVGG->AddEntry(&dummy_vgg_unpol, "VGG unpolarized", "l");
    } else if (mode == XSecPanelMode::PosOnly) {
        legData->AddEntry(&dummy_pos, "data + helicity", "lep");
        legKM->AddEntry(&dummy_bh,     "BH",            "l");
        legKM->AddEntry(&dummy_km_pos, "KM + helicity", "l");
        legVGG->AddEntry(&dummy_vgg_pos, "VGG + helicity", "l");
    } else if (mode == XSecPanelMode::NegOnly) {
        legData->AddEntry(&dummy_neg, "data - helicity", "lep");
        legKM->AddEntry(&dummy_bh,     "BH",            "l");
        legKM->AddEntry(&dummy_km_neg, "KM - helicity", "l");
        legVGG->AddEntry(&dummy_vgg_neg, "VGG - helicity", "l");
    }

    legData->Draw();
    legKM->Draw();
    legVGG->Draw();

    for (int r = 0; r < nrows; ++r) {
        const Range &t_range = t_slice[r];

        for (int cc = 0; cc < ncols; ++cc) {
            const Range &q2_range = q2_slice[cc];

            pGrid->cd(r * ncols + cc + 1);
            gPad->SetGrid(1, 1);
            gPad->SetTopMargin(0.12);
            gPad->SetBottomMargin(0.18);
            gPad->SetLeftMargin(0.16);
            gPad->SetRightMargin(0.10);
            gPad->SetLogy(true);

            QTKey key(q2_range, t_range);
            auto it_bin = bins_for_xB.find(key);

            const BinData *bin_ptr = nullptr;

            if (it_bin != bins_for_xB.end()) {
                bin_ptr = &(it_bin->second);
            }

            double ymin_canvas = 1e-4;
            double ymax_canvas = 1.0;

            {
                auto yr = compute_yrange_for_bin(bin_ptr, theory, mode);
                ymin_canvas = yr.first;
                ymax_canvas = yr.second;
            }

            TH1 *frame = gPad->DrawFrame(0.0, ymin_canvas, 360.0, ymax_canvas);
            frame->GetXaxis()->SetTitle("#phi (deg)");
            frame->GetYaxis()->SetTitle("d^{4}#sigma / (dx_{B} dQ^{2} d|t| d#phi)");
            frame->GetXaxis()->CenterTitle();
            frame->GetYaxis()->CenterTitle();
            frame->GetXaxis()->SetNdivisions(505);
            frame->GetXaxis()->SetTitleSize(0.060);
            frame->GetYaxis()->SetTitleSize(0.060);
            frame->GetXaxis()->SetLabelSize(0.048);
            frame->GetYaxis()->SetLabelSize(0.048);
            frame->GetXaxis()->SetTitleOffset(1.10);
            frame->GetYaxis()->SetTitleOffset(1.35);

            TLatex lab;
            lab.SetNDC();
            lab.SetTextSize(cellLabelSize);
            lab.SetTextAlign(11);
            lab.SetTextFont(42);
            lab.DrawLatex(
                0.14, 0.93,
                Form("Q^{2} in (%.2f, %.2f), |t| in (%.2f, %.2f)",
                     q2_range.first, q2_range.second,
                     t_range.first,  t_range.second)
            );

            if (!bin_ptr) continue;

            BinData bin = *bin_ptr;

            auto sort_by_phi = [](std::vector<Point> &v) {
                std::sort(v.begin(), v.end(),
                          [](const Point &a, const Point &b) {
                              return a.phi < b.phi;
                          });
            };

            sort_by_phi(bin.unpol);
            sort_by_phi(bin.pos);
            sort_by_phi(bin.neg);

            if (mode == XSecPanelMode::All || mode == XSecPanelMode::UnpolOnly) {
                TGraphErrors *g_unpol = make_xsec_graph(bin.unpol, 20, kBlack);

                if (g_unpol) {
                    g_unpol->Draw("P SAME");
                }
            }

            if (mode == XSecPanelMode::All || mode == XSecPanelMode::PosOnly) {
                TGraphErrors *g_pos = make_xsec_graph(bin.pos, 24, kRed + 1);

                if (g_pos) {
                    g_pos->Draw("P SAME");
                }
            }

            if (mode == XSecPanelMode::All || mode == XSecPanelMode::NegOnly) {
                TGraphErrors *g_neg = make_xsec_graph(bin.neg, 25, kBlue + 1);

                if (g_neg) {
                    g_neg->Draw("P SAME");
                }
            }

            if (bin.have_theory_row) {
                auto it_th = theory.find(bin.theory_row);

                if (it_th != theory.end()) {
                    const TheoryCurves &tc = it_th->second;

                    auto draw_curve = [&](const std::vector<double> &ys,
                                          int lstyle,
                                          int lcolor) {
                        if (tc.phi_deg.size() != ys.size() || ys.empty()) return;

                        int M = (int)ys.size();

                        std::vector<double> xp(M);
                        std::vector<double> yp(M);

                        for (int i = 0; i < M; ++i) {
                            xp[i] = tc.phi_deg[i];
                            yp[i] = ys[i];
                        }

                        TGraph *gth = new TGraph(M, xp.data(), yp.data());
                        gth->SetLineStyle(lstyle);
                        gth->SetLineWidth(2);
                        gth->SetLineColor(lcolor);
                        gth->Draw("L SAME");
                    };

                    if (mode == XSecPanelMode::All) {
                        draw_curve(tc.bh_unpol, 2, kGreen + 2);
                        draw_curve(tc.km_unpol, 1, kMagenta + 1);
                        draw_curve(tc.km_pos,   2, kMagenta + 1);
                        draw_curve(tc.km_neg,   3, kMagenta + 1);
                        draw_curve(tc.vgg_unpol, 1, kOrange + 7);
                        draw_curve(tc.vgg_pos,   2, kOrange + 7);
                        draw_curve(tc.vgg_neg,   3, kOrange + 7);
                    } else if (mode == XSecPanelMode::UnpolOnly) {
                        draw_curve(tc.bh_unpol, 2, kGreen + 2);
                        draw_curve(tc.km_unpol, 1, kMagenta + 1);
                        draw_curve(tc.vgg_unpol, 1, kOrange + 7);
                    } else if (mode == XSecPanelMode::PosOnly) {
                        draw_curve(tc.bh_unpol, 2, kGreen + 2);
                        draw_curve(tc.km_pos,   2, kMagenta + 1);
                        draw_curve(tc.vgg_pos,  2, kOrange + 7);
                    } else if (mode == XSecPanelMode::NegOnly) {
                        draw_curve(tc.bh_unpol, 2, kGreen + 2);
                        draw_curve(tc.km_neg,   3, kMagenta + 1);
                        draw_curve(tc.vgg_neg,  3, kOrange + 7);
                    }
                }
            }
        }
    }

    std::ostringstream fname;
    fname << "cross_sections_";

    if (mode == XSecPanelMode::UnpolOnly) {
        fname << "unpol_";
    } else if (mode == XSecPanelMode::PosOnly) {
        fname << "pos_";
    } else if (mode == XSecPanelMode::NegOnly) {
        fname << "neg_";
    }

    fname << canonical_period_dir(label)
          << "_xB_" << xb_idx_for_name << ".png";

    fs::path outpath = outdir / fname.str();
    c->SaveAs(outpath.string().c_str());

    delete c;
}

// -----------------------------------------------------------------------------
// Main plotting function
// -----------------------------------------------------------------------------

bool plot_cross_sections_for_label(const std::string &csv_main,
                                   const std::string &label,
                                   const std::string &theory_json_root,
                                   const std::string &out_root_dir) {
    std::ifstream ifs(csv_main);

    if (!ifs) {
        std::cerr << "[cross_sections] ERROR: cannot open " << csv_main
                  << " for plotting.\n";
        return false;
    }

    std::vector<std::string> lines;
    std::string line;

    while (std::getline(ifs, line)) {
        lines.push_back(line);
    }

    ifs.close();

    if (lines.empty()) {
        std::cerr << "[cross_sections] ERROR: CSV " << csv_main
                  << " is empty in plot_cross_sections_for_label.\n";
        return false;
    }

    std::vector<std::string> header = split_csv_line(lines[0]);

    int c_xb_min = -1;
    int c_xb_max = -1;
    int c_q2_min = -1;
    int c_q2_max = -1;
    int c_t_min  = -1;
    int c_t_max  = -1;

    try {
        c_xb_min = find_col(header, "xBmin");
        c_xb_max = find_col(header, "xBmax");
        c_q2_min = find_col(header, "Q2min");
        c_q2_max = find_col(header, "Q2max");
        c_t_min  = find_col(header, "t_abs_min");
        c_t_max  = find_col(header, "t_abs_max");
    } catch (const std::exception &e) {
        std::cerr << "[cross_sections] FATAL: " << e.what()
                  << " in plot_cross_sections_for_label.\n";
        return false;
    }

    int c_xb_idx = find_col_optional(header, "xB index");

    int c_phiavg = find_col_optional(header, "phiavg, " + label);
    int c_phimin = -1;
    int c_phimax = -1;

    if (c_phiavg < 0) {
        c_phimin = find_col_optional(header, "phimin");
        c_phimax = find_col_optional(header, "phimax");

        if (c_phimin < 0 || c_phimax < 0) {
            std::cerr << "[cross_sections] FATAL: no phiavg or phimin/phimax "
                      << "available for label " << label << ".\n";
            return false;
        }
    }

    int c_xs_unpol = find_col_optional(
        header, "cross sections, ep->epg, exp, " + label + ", unpol");

    int c_xs_pos = find_col_optional(
        header, "cross sections, ep->epg, exp, " + label + ", pos");

    int c_xs_neg = find_col_optional(
        header, "cross sections, ep->epg, exp, " + label + ", neg");

    if (c_xs_unpol < 0 && c_xs_pos < 0 && c_xs_neg < 0) {
        std::cerr << "[cross_sections] INFO: no cross section columns for label "
                  << label << "; nothing to plot.\n";
        return true;
    }

    if (c_xs_unpol < 0) {
        std::cerr << "[cross_sections] FATAL: missing unpolarized cross section column "
                  << "for label " << label << ".\n";
        return false;
    }

    const bool has_hel = has_helicity_resolved_cross_sections(label);

    if (has_hel && (c_xs_pos < 0 || c_xs_neg < 0)) {
        std::cerr << "[cross_sections] FATAL: incomplete helicity-resolved cross section columns "
                  << "for label " << label << " (need unpol/pos/neg).\n";
        return false;
    }

    std::map<Range, XSGroupByXB> by_xb;

    for (size_t row = 1; row < lines.size(); ++row) {
        if (lines[row].empty()) continue;

        std::vector<std::string> fields = split_csv_line(lines[row]);

        if (fields.size() != header.size()) continue;

        double xbmin = std::atof(trim(unquote(fields[c_xb_min])).c_str());
        double xbmax = std::atof(trim(unquote(fields[c_xb_max])).c_str());
        double q2min = std::atof(trim(unquote(fields[c_q2_min])).c_str());
        double q2max = std::atof(trim(unquote(fields[c_q2_max])).c_str());
        double tmin  = std::atof(trim(unquote(fields[c_t_min])).c_str());
        double tmax  = std::atof(trim(unquote(fields[c_t_max])).c_str());

        Range xb_range(xbmin, xbmax);
        Range q2_range(q2min, q2max);
        Range t_range(tmin,  tmax);

        double phi = 0.0;

        if (c_phiavg >= 0) {
            phi = std::atof(trim(unquote(fields[c_phiavg])).c_str());
        } else {
            double pmin = std::atof(trim(unquote(fields[c_phimin])).c_str());
            double pmax = std::atof(trim(unquote(fields[c_phimax])).c_str());
            phi = 0.5 * (pmin + pmax);
        }

        Triple xs_unpol = parse_tuple3(fields[c_xs_unpol]);
        Triple xs_pos{0.0, 0.0, 0.0};
        Triple xs_neg{0.0, 0.0, 0.0};

        if (has_hel) {
            xs_pos = parse_tuple3(fields[c_xs_pos]);
            xs_neg = parse_tuple3(fields[c_xs_neg]);
        }

        if (xs_unpol.value <= 0.0 &&
            (!has_hel || (xs_pos.value <= 0.0 && xs_neg.value <= 0.0))) {
            continue;
        }

        XSGroupByXB &group = by_xb[xb_range];

        if (group.xb_index < 0 && c_xb_idx >= 0) {
            group.xb_index = std::atoi(trim(unquote(fields[c_xb_idx])).c_str());
        }

        QTKey key(q2_range, t_range);
        BinData &bin = group.bins[key];

        if (!bin.have_theory_row) {
            bin.theory_row      = row;
            bin.have_theory_row = true;
        }

        auto add_point = [&](const Triple &xs, std::vector<Point> &vec) {
            if (xs.value <= 0.0) return;

            Point p;
            p.phi    = phi;
            p.xs     = xs.value;
            p.xs_err = xs.stat;
            vec.push_back(p);
        };

        add_point(xs_unpol, bin.unpol);

        if (has_hel) {
            add_point(xs_pos, bin.pos);
            add_point(xs_neg, bin.neg);
        }
    }

    if (by_xb.empty()) {
        std::cerr << "[cross_sections] WARNING: no xsec data found for label "
                  << label << " to plot.\n";
        return true;
    }

    std::map<size_t, TheoryCurves> theory = load_theory_for_label(label, theory_json_root);

    gROOT->SetBatch(true);
    gStyle->SetOptStat(0);
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetLabelFont(42, "XYZ");
    gStyle->SetTextFont(42);
    gStyle->SetLineWidth(2);
    gStyle->SetFrameLineWidth(2);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    fs::path outdir = fs::path(out_root_dir) / canonical_period_dir(label);
    ensure_dir(outdir);

    int xb_canvas_counter = 0;

    for (const auto &kv_xb : by_xb) {
        const Range &xb_range = kv_xb.first;
        const XSGroupByXB &group = kv_xb.second;
        const auto &bins_for_xB = group.bins;

        if (bins_for_xB.empty()) continue;

        std::set<Range> q2_set;
        std::set<Range> t_set;

        for (const auto &kv : bins_for_xB) {
            const QTKey &qt = kv.first;
            q2_set.insert(qt.first);
            t_set.insert(qt.second);
        }

        std::vector<Range> q2_slice(q2_set.begin(), q2_set.end());
        std::vector<Range> t_slice(t_set.begin(), t_set.end());

        if (q2_slice.empty() || t_slice.empty()) continue;

        int ncols = (int)q2_slice.size();
        int nrows = (int)t_slice.size();
        int nPads = ncols * nrows;

        int xb_idx_for_name =
            (group.xb_index >= 0 ? group.xb_index : xb_canvas_counter);

        if (has_hel) {
            make_xsec_canvas_for_mode(label, xb_range, group,
                                      q2_slice, t_slice,
                                      theory, outdir,
                                      xb_idx_for_name,
                                      XSecPanelMode::All,
                                      ncols, nrows, nPads);
        }

        make_xsec_canvas_for_mode(label, xb_range, group,
                                  q2_slice, t_slice,
                                  theory, outdir,
                                  xb_idx_for_name,
                                  XSecPanelMode::UnpolOnly,
                                  ncols, nrows, nPads);

        if (has_hel) {
            make_xsec_canvas_for_mode(label, xb_range, group,
                                      q2_slice, t_slice,
                                      theory, outdir,
                                      xb_idx_for_name,
                                      XSecPanelMode::PosOnly,
                                      ncols, nrows, nPads);

            make_xsec_canvas_for_mode(label, xb_range, group,
                                      q2_slice, t_slice,
                                      theory, outdir,
                                      xb_idx_for_name,
                                      XSecPanelMode::NegOnly,
                                      ncols, nrows, nPads);
        }

        ++xb_canvas_counter;
    }

    std::cout << "[cross_sections] Plotted cross sections for label "
              << label << " into " << outdir.string() << "\n";

    return true;
}