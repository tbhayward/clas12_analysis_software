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
//   - It reads correction factors already stored in the CSV:
//       Frad, <energy>
//       Fbin, <energy>
//       bin_volume, <energy>
//
//   - It computes:
//       sigma = Y_unfolded * Frad * Fbin / (L * bin_volume)
//
//   - It writes:
//       cross sections, ep->epg, exp, <label>, <helicity>
//
//   - It also fills integrated luminosity columns in the CSV. For combined
//     labels, luminosity remains row-dependent and includes only member periods
//     whose acceptance is nonzero in that row. This avoids normalizing a combined
//     yield by charge from a period that did not contribute acceptance in the bin.
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
    if (!fs::exists(p)) {
        fs::create_directories(p);
    }
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
    Triple out;
    out.value = 0.0;
    out.stat  = 0.0;
    out.sys   = 0.0;

    std::string s = strip_all_outer_quotes(cell);
    s = trim(s);

    if (s.empty()) {
        return out;
    }

    if (s.front() == '(' && s.back() == ')') {
        s = s.substr(1, s.size() - 2);
        s = trim(s);

        if (s.empty()) {
            return out;
        }
    }

    std::vector<std::string> parts;
    std::string token;

    for (char c : s) {
        if (c == ',') {
            parts.push_back(trim(token));
            token.clear();
        } else {
            token.push_back(c);
        }
    }

    if (!token.empty()) {
        parts.push_back(trim(token));
    }

    auto to_double_or_zero = [](const std::string &str) -> double {
        if (str.empty()) {
            return 0.0;
        }

        return std::atof(str.c_str());
    };

    if (!parts.empty()) {
        out.value = to_double_or_zero(parts[0]);
    }

    if (parts.size() > 1U) {
        out.stat = to_double_or_zero(parts[1]);
    }

    if (parts.size() > 2U) {
        out.sys = to_double_or_zero(parts[2]);
    }

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

// -----------------------------------------------------------------------------
// Luminosity helpers
// -----------------------------------------------------------------------------

static Triple load_rga_lumi_file(const std::string &path,
                                 bool unpolarized_total_from_pos_plus_neg,
                                 bool use_scaled_columns_3_to_5_for_unpolarized,
                                 double columns_3_to_5_charge_sum_scale) {
    std::ifstream ifs(path);

    if (!ifs) {
        std::cerr << "[cross_sections] FATAL: cannot open lumi file: "
                  << path << "\n";
        throw std::runtime_error("cannot open lumi file");
    }

    double sum_total_col = 0.0;
    double sum_pos       = 0.0;
    double sum_neg       = 0.0;
    double sum_col5      = 0.0;

    std::string line;
    size_t n_lines = 0;

    while (std::getline(ifs, line)) {
        std::string s = trim(line);

        if (s.empty()) continue;
        if (!s.empty() && s[0] == '#') continue;

        std::vector<std::string> fields;
        std::string field;
        std::istringstream iss(s);

        while (std::getline(iss, field, ',')) {
            fields.push_back(trim(field));
        }

        if (fields.size() < 4) {
            std::cerr << "[cross_sections] WARNING: lumi file " << path
                      << " has a line with fewer than 4 columns, skipping: "
                      << s << "\n";
            continue;
        }

        if (use_scaled_columns_3_to_5_for_unpolarized && fields.size() < 5) {
            std::cerr << "[cross_sections] WARNING: lumi file " << path
                      << " has a line with fewer than 5 columns in scaled columns-3-to-5 mode, skipping: "
                      << s << "\n";
            continue;
        }

        const double total_col = std::atof(fields[1].c_str());
        const double pos       = std::atof(fields[2].c_str());
        const double neg       = std::atof(fields[3].c_str());
        const double col5      = (fields.size() >= 5) ? std::atof(fields[4].c_str()) : 0.0;

        sum_total_col += total_col;
        sum_pos       += pos;
        sum_neg       += neg;
        sum_col5      += col5;
        ++n_lines;
    }

    double final_total = sum_total_col;
    std::string source = "column 2";

    if (use_scaled_columns_3_to_5_for_unpolarized) {
        final_total = columns_3_to_5_charge_sum_scale * (sum_pos + sum_neg + sum_col5);
        source = "scaled columns 3+4+5";
    } else if (unpolarized_total_from_pos_plus_neg) {
        final_total = sum_pos + sum_neg;
        source = "pos+neg columns";
    }

    std::cout << "[cross_sections] Loaded lumi from " << path
              << " over " << n_lines << " runs: "
              << "unpolarized_total=" << final_total
              << " total_col=" << sum_total_col
              << " pos=" << sum_pos
              << " neg=" << sum_neg
              << " col5=" << sum_col5
              << "  [unpolarized source: " << source;

    if (use_scaled_columns_3_to_5_for_unpolarized) {
        std::cout << ", scale=" << columns_3_to_5_charge_sum_scale;
    }

    std::cout << "]\n";

    Triple out;
    out.value = final_total;
    out.stat  = sum_pos;
    out.sys   = sum_neg;
    return out;
}


LumiMap build_lumi_map() {
    LumiBuildOptions options;
    options.use_second_column_charge_for_all_unpolarized = true;
    return build_lumi_map(options);
}

LumiMap build_lumi_map(const LumiBuildOptions &options) {
    LumiMap m;
    const std::string base = "imports/integrated_luminosity";

    const bool use_col2_for_all =
        options.use_second_column_charge_for_all_unpolarized;

    const bool use_scaled_fa18_sp19 =
        options.use_columns_3_to_5_charge_sum_scaled_for_fa18_sp19_unpolarized;

    std::cout << "[cross_sections] build_lumi_map charge convention: ";

    if (use_scaled_fa18_sp19) {
        std::cout << "Fa18/Sp19 use scale * (columns 3+4+5), scale="
                  << options.columns_3_to_5_charge_sum_scale
                  << "; Sp18 uses column 2.\n";
    } else if (use_col2_for_all) {
        std::cout << "using column 2 for unpolarized accumulated charge for all periods.\n";
    } else {
        std::cout << "legacy mixed mode: Sp18 uses column 2; Fa18 and Sp19 use columns 3+4.\n";
    }

    // These settings control only Triple.value, i.e. the unpolarized luminosity.
    // Triple.stat and Triple.sys are always filled from columns 3 and 4,
    // respectively, so polarized cross sections remain normalized with the
    // helicity-specific accumulated charges.
    // Spring 2018 is always forced to column 2 for unpolarized normalization.
    const bool fa18_unpol_from_pos_neg = (!use_col2_for_all && !use_scaled_fa18_sp19);
    const bool sp19_unpol_from_pos_neg = (!use_col2_for_all && !use_scaled_fa18_sp19);
    const bool sp18_unpol_from_pos_neg = false;

    const bool fa18_use_scaled_cols_3_to_5 = use_scaled_fa18_sp19;
    const bool sp19_use_scaled_cols_3_to_5 = use_scaled_fa18_sp19;
    const bool sp18_use_scaled_cols_3_to_5 = false;

    try {
        m["Fa18 Inb"]      = load_rga_lumi_file(base + "/rga_fa18_inb.txt",
                                                fa18_unpol_from_pos_neg,
                                                fa18_use_scaled_cols_3_to_5,
                                                options.columns_3_to_5_charge_sum_scale);

        m["Fa18 Out"]      = load_rga_lumi_file(base + "/rga_fa18_out.txt",
                                                fa18_unpol_from_pos_neg,
                                                fa18_use_scaled_cols_3_to_5,
                                                options.columns_3_to_5_charge_sum_scale);

        m["Fa18 Inb Supp"] = Triple{0.0, 0.0, 0.0};

        m["Sp18 Inb"]      = load_rga_lumi_file(base + "/rga_sp18_inb.txt",
                                                sp18_unpol_from_pos_neg,
                                                sp18_use_scaled_cols_3_to_5,
                                                options.columns_3_to_5_charge_sum_scale);

        m["Sp18 Out"]      = load_rga_lumi_file(base + "/rga_sp18_out.txt",
                                                sp18_unpol_from_pos_neg,
                                                sp18_use_scaled_cols_3_to_5,
                                                options.columns_3_to_5_charge_sum_scale);

        m["Sp19 Inb"]      = load_rga_lumi_file(base + "/rga_sp19_inb.txt",
                                                sp19_unpol_from_pos_neg,
                                                sp19_use_scaled_cols_3_to_5,
                                                options.columns_3_to_5_charge_sum_scale);
    } catch (const std::exception &e) {
        std::cerr << "[cross_sections] FATAL in build_lumi_map: "
                  << e.what() << "\n";
        throw;
    }

    auto sum_labels = [&](const std::vector<std::string> &keys) -> Triple {
        Triple r{0.0, 0.0, 0.0};

        for (const auto &k : keys) {
            auto it = m.find(k);

            if (it == m.end()) {
                std::cerr << "[cross_sections] ERROR: missing lumi for \""
                          << k << "\" while building combined groups.\n";
                continue;
            }

            r.value += it->second.value;
            r.stat  += it->second.stat;
            r.sys   += it->second.sys;
        }

        return r;
    };

    // Supplemental is intentionally not included in the combined groups used for
    // the production cross-section normalization.
    m["Fa18"]     = sum_labels({"Fa18 Inb", "Fa18 Out"});
    m["Sp18"]     = sum_labels({"Sp18 Inb", "Sp18 Out"});
    m["10.6 GeV"] = sum_labels({"Fa18 Inb", "Fa18 Out", "Sp18 Inb", "Sp18 Out"});
    m["10.2 GeV"] = sum_labels({"Sp19 Inb"});

    std::cout << "[cross_sections] build_lumi_map summary:\n"
              << "  Fa18 Inb  value=" << m["Fa18 Inb"].value
              << " pos=" << m["Fa18 Inb"].stat
              << " neg=" << m["Fa18 Inb"].sys << "\n"
              << "  Fa18 Out  value=" << m["Fa18 Out"].value
              << " pos=" << m["Fa18 Out"].stat
              << " neg=" << m["Fa18 Out"].sys << "\n"
              << "  Sp18 Inb  value=" << m["Sp18 Inb"].value
              << " pos=" << m["Sp18 Inb"].stat
              << " neg=" << m["Sp18 Inb"].sys << "\n"
              << "  Sp18 Out  value=" << m["Sp18 Out"].value
              << " pos=" << m["Sp18 Out"].stat
              << " neg=" << m["Sp18 Out"].sys << "\n"
              << "  Sp19 Inb  value=" << m["Sp19 Inb"].value
              << " pos=" << m["Sp19 Inb"].stat
              << " neg=" << m["Sp19 Inb"].sys << "\n"
              << "  Fa18      value=" << m["Fa18"].value
              << " pos=" << m["Fa18"].stat
              << " neg=" << m["Fa18"].sys << "\n"
              << "  Sp18      value=" << m["Sp18"].value
              << " pos=" << m["Sp18"].stat
              << " neg=" << m["Sp18"].sys << "\n"
              << "  10.6 GeV  value=" << m["10.6 GeV"].value
              << " pos=" << m["10.6 GeV"].stat
              << " neg=" << m["10.6 GeV"].sys << "\n"
              << "  10.2 GeV  value=" << m["10.2 GeV"].value
              << " pos=" << m["10.2 GeV"].stat
              << " neg=" << m["10.2 GeV"].sys << "\n";

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

    auto period_has_acceptance = [&](const std::string &period,
                                     const std::vector<std::string> &fields) -> bool {
        auto it = acc_col_idx.find(period);

        if (it == acc_col_idx.end()) return false;

        const Triple a = parse_tuple3(fields[it->second]);
        return a.value > 0.0;
    };

    auto lumi_for_label_row = [&](const std::string &L,
                                  const std::vector<std::string> &fields) -> Triple {
        if (!is_combined_label(L)) {
            auto it = lumi_map.find(L);

            if (it == lumi_map.end()) return Triple{0.0, 0.0, 0.0};

            return it->second;
        }

        Triple out{0.0, 0.0, 0.0};

        for (const auto &m : combined_members_for(L)) {
            if (!period_has_acceptance(m, fields)) continue;

            auto itLm = lumi_map.find(m);

            if (itLm == lumi_map.end()) continue;

            out.value += itLm->second.value;
            out.stat  += itLm->second.stat;
            out.sys   += itLm->second.sys;
        }

        return out;
    };

    auto write_lumi_columns_for_row = [&](std::vector<std::string> &fields) {
        for (const auto &L : labels) {
            const std::string name = lumi_col_for_label(L);

            if (name.empty()) continue;

            auto itc = lumi_col_idx.find(name);

            if (itc == lumi_col_idx.end()) continue;

            Triple lum = lumi_for_label_row(L, fields);
            fields[itc->second] = tuple3_to_cell(lum.value, lum.stat, lum.sys);
        }
    };

    std::vector<std::string> out_lines;
    out_lines.reserve(lines.size());
    out_lines.push_back(lines[0]);

    std::cout << "[cross_sections] compute_cross_sections: data rows = "
              << n_data_rows << "\n";

    std::cout << "[cross_sections] NOTE: no imports/efficiency.json correction is applied here. "
              << "Current-efficiency and eppi0 normalization corrections are already upstream.\n";

    std::cout << "[cross_sections] NOTE: combined-label luminosities are row-dependent and "
              << "gated by nonzero member-period acceptance.\n";

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

        write_lumi_columns_for_row(fields);

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

        for (const auto &L : labels) {
            const Triple Lumi = lumi_for_label_row(L, fields);
            const bool use_10p2 = (L == "Sp19 Inb" || L == "10.2 GeV");

            const Triple &Frad = frad;
            const Triple &Fbin = fbin;

            if (Frad.value <= 0.0 || Fbin.value <= 0.0) continue;

            for (const auto &h : cross_section_helicities_for_label(L)) {
                auto it = colmap.find(L + "|" + h);

                if (it == colmap.end()) continue;

                const int iy = it->second.yield_idx;
                const int ix = it->second.xs_idx;

                const Triple Y = parse_tuple3(fields[iy]);

                if (Y.value <= 0.0) continue;

                const int c_vbin = use_10p2 ? c_vbin_102 : c_vbin_106;
                const int c_cubic_vbin = use_10p2 ? c_cubic_vbin_102 : c_cubic_vbin_106;
                const std::string energy_tag = use_10p2 ? "10.2 GeV" : "10.6 GeV";

                const Triple Vbin = parse_required_triple_cell(
                    fields,
                    c_vbin,
                    "bin_volume, " + energy_tag,
                    "pass-2 CSV row " + std::to_string(row) + ", label " + L + ", helicity " + h
                );

                // Validate the diagnostic cubic volume for the same energy. It is not
                // used in sigma, but it confirms this row came from the new schema
                // and lets the CSV expose allowed/cubic volume ratios for debugging.
                (void)parse_required_triple_cell(
                    fields,
                    c_cubic_vbin,
                    "cubic bin_volume, " + energy_tag,
                    "pass-2 CSV row " + std::to_string(row) + ", label " + L + ", helicity " + h
                );

                double lumi_val = 0.0;

                if (h == "unpol") {
                    lumi_val = Lumi.value;
                } else if (h == "pos") {
                    lumi_val = Lumi.stat;
                } else if (h == "neg") {
                    lumi_val = Lumi.sys;
                }

                if (lumi_val <= 0.0) continue;

                const double denom = lumi_val * Vbin.value;

                if (denom <= 0.0) continue;

                const double sigma = Y.value * Frad.value * Fbin.value / denom;

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