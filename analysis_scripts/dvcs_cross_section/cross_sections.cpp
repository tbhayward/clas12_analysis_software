// cross_sections.cpp
// -----------------------------------------------------------------------------
// Cross section calculator and plotting with BH / KM / VGG theory curves.
//
// This module:
//   1) Updates dvcs_pass2_analysis.csv:
//        - Fills integrated luminosity columns for all rows (where such
//          columns already exist).
//        - Computes cross sections for unpol/pos/neg per label where both
//          the acceptance corrected yield column and the cross section
//          column exist.
//        - Stores them as "(value, stat, sys)" with sys = 0.
//        - Propagates stat uncertainties from the acceptance corrected yield,
//          Frad, and (optionally) Fbin.
//        - By default, uses Lee's original Frad/Fbin columns ("Frad", "Fbin")
//          (copied from all_bin_v3.csv) but **inverts** those values internally
//          to recover the standard multiplicative Frad and bin-centering
//          corrections. The **bin volume** is always taken from the pass-2
//          per-energy bin_volume columns ("bin_volume, 10.6 GeV" /
//          "bin_volume, 10.2 GeV").
//
//          Summary of behavior:
//            * kUseLeeFradFbin == true:
//                 - Use Lee's "Frad" and "Fbin" (inverted) as multiplicative
//                   Frad and Fbin bin-centering corrections.
//                 - Use bin_volume from pass-2 per energy (10.6 / 10.2).
//            * kUseLeeFradFbin == false:
//                 - Use pass-2 per-energy Frad and Fbin columns:
//                      "Frad, 10.6 GeV" / "Frad, 10.2 GeV",
//                      "Fbin, 10.6 GeV" / "Fbin, 10.2 GeV".
//                 - Use bin_volume from pass-2 per energy.
//
//   2) Assumes the BH / KM / VGG theory curves as a function of phi have
//      already been computed at 10.6 GeV and 10.2 GeV and stored in JSON files:
//          output/jsons/cross_sections/10.6_GeV/xs_phi_all.json
//          output/jsons/cross_sections/10.2_GeV/xs_phi_all.json
//      When plotting *any* label (Fa18 Inb, Fa18 Out, Sp18 Inb, Sp18 Out,
//      Sp19 Inb, Fa18, Sp18, 10.6 GeV, 10.2 GeV, ...), we select the theory
//      JSON based on the beam energy (10.6 vs 10.2).
//   3) For each label, plots xB canvases where:
//        - One canvas per xB bin.
//        - Each canvas is laid out as N_t rows (|t| bins) by N_Q2 columns
//          (Q^2 bins).
//        - Each subpad (a fixed Q^2 and |t| bin) shows:
//            * unpolarized, + helicity, and - helicity cross sections vs phi
//              (with errors) from the CSV.
//            * BH (unpol), KM (unpol/+/−), VGG (unpol/+/−) curves vs phi
//              from the theory JSON (10.6 or 10.2 GeV as appropriate).
//        - Y axis is in log scale, but each subpad has its own dynamic
//          y-range determined from its data + theory.
//        - A single title plus three separate legend boxes (Data / BH+KM / VGG)
//          are drawn in the top pad of the canvas.
//        - Empty (Q^2, |t|) bins show only an annotated frame.
//   4) This file does not add new columns to the CSV; it only updates
//      existing ones. If required columns are missing, it fails with a clear
//      error; optional combined-group lumi columns such as
//      "integrated luminosity, 10.2 GeV (nC)" are treated as optional and
//      simply not filled if absent.
//   5) In addition to the full "all helicities" canvases, we also make
//      three extra versions for each xB bin:
//        - Unpolarized only (data unpol + BH/KM/VGG unpol).
//        - + helicity only (data + helicity + BH unpol + KM/VGG +).
//        - - helicity only (data - helicity + BH unpol + KM/VGG -).
//      These are saved with suffixes _unpol, _pos, and _neg in the filename.
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

// ROOT includes (for plotting part only)
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
using json = nlohmann::json;

// Convenience typedef for (min,max) ranges
using Range = std::pair<double,double>;

// -----------------------------------------------------------------------------
// Configuration toggles
// -----------------------------------------------------------------------------

static const bool kUseLeeFradFbin = true;

// -----------------------------------------------------------------------------
// Lee pass-1 Frad/Fbin table loader
// -----------------------------------------------------------------------------

struct LeeFradFbinRow {
    double frad;
    double fbin;
};

using LeeFradFbinTable = std::vector<LeeFradFbinRow>;

// Path to Lee's all_bin_v3.csv (relative to your run directory).
static const std::string kLeeAllBinCsvPath = "imports/all_bin_v3.csv";

// -----------------------------------------------------------------------------
// Forward declarations for CSV helpers used before their full definitions
// -----------------------------------------------------------------------------

static std::vector<std::string> split_csv_line(const std::string &line);
static std::string trim(const std::string &s);
static std::string unquote(const std::string &s);
static int find_col(const std::vector<std::string> &header,
                    const std::string &target);

static LeeFradFbinTable load_lee_frad_fbin_table(const std::string &path) {
    std::ifstream ifs(path);
    if (!ifs) {
        std::cerr << "[cross_sections] FATAL: cannot open Lee Frad/Fbin file: "
                  << path << "\n";
        throw std::runtime_error("cannot open Lee Frad/Fbin file: " + path);
    }

    std::vector<std::string> lines;
    std::string line;
    while (std::getline(ifs, line)) {
        lines.push_back(line);
    }
    ifs.close();

    if (lines.empty()) {
        std::cerr << "[cross_sections] FATAL: Lee Frad/Fbin file "
                  << path << " is empty.\n";
        throw std::runtime_error("Lee Frad/Fbin file empty: " + path);
    }

    std::vector<std::string> header = split_csv_line(lines[0]);

    int c_frad = -1;
    int c_fbin = -1;
    try {
        c_frad = find_col(header, "Frad");
        c_fbin = find_col(header, "Fbin");
    } catch (const std::exception &e) {
        std::cerr << "[cross_sections] FATAL while reading Lee Frad/Fbin file "
                  << path << ": " << e.what() << "\n";
        throw;
    }

    LeeFradFbinTable table;
    table.reserve(lines.size() - 1);

    for (size_t row = 1; row < lines.size(); ++row) {
        if (lines[row].empty()) {
            table.push_back({0.0, 0.0});
            continue;
        }

        std::vector<std::string> fields = split_csv_line(lines[row]);
        if (fields.size() != header.size()) {
            std::ostringstream oss;
            oss << "Line " << row << " in " << path
                << " has " << fields.size()
                << " columns, expected " << header.size();
            std::cerr << "[cross_sections] FATAL: " << oss.str() << "\n";
            throw std::runtime_error(oss.str());
        }

        std::string frad_str = trim(unquote(fields[c_frad]));
        std::string fbin_str = trim(unquote(fields[c_fbin]));

        double frad_val = 0.0;
        double fbin_val = 0.0;

        if (!frad_str.empty()) {
            frad_val = std::atof(frad_str.c_str());
        }
        if (!fbin_str.empty()) {
            fbin_val = std::atof(fbin_str.c_str());
        }

        table.push_back({frad_val, fbin_val});
    }

    std::cout << "[cross_sections] Loaded Lee Frad/Fbin table from "
              << path << " with " << table.size() << " data rows.\n";

    return table;
}

// -----------------------------------------------------------------------------
// Period directory helper
// -----------------------------------------------------------------------------

static std::string canonical_period_dir(const std::string &label) {
    if (label == "Fa18 Inb")       return "Fa18_Inb";
    if (label == "Fa18 Out")       return "Fa18_Out";
    if (label == "Fa18 Inb Supp")  return "Fa18_Inb_Supp";
    if (label == "Sp18 Inb")       return "Sp18_Inb";
    if (label == "Sp18 Out")       return "Sp18_Out";
    if (label == "Sp19 Inb")       return "Sp19_Inb";
    if (label == "Fa18")           return "Fa18";
    if (label == "Sp18")           return "Sp18";
    if (label == "10.6 GeV")       return "10.6_GeV";
    if (label == "10.2 GeV")       return "10.2_GeV";

    std::string out = label;
    std::replace(out.begin(), out.end(), ' ', '_');
    return out;
}

static std::string yield_label_for(const std::string &label) {
    if (label == "10.6 GeV") {
        return "2018 (10.6 GeV)";
    }
    if (label == "10.2 GeV") {
        return "Sp19 Inb";
    }
    return label;
}

static double beam_energy_for_label(const std::string &label) {
    if (label == "Sp19 Inb" || label == "10.2 GeV") {
        return 10.2;
    }
    return 10.6;
}

// -----------------------------------------------------------------------------
// CSV helpers
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

// Triple is defined in cross_sections.h:
//   struct Triple { double value; double stat; double sys; };

static Triple parse_tuple3(const std::string &cell) {
    std::string s = trim(unquote(cell));
    Triple out{0.0, 0.0, 0.0};
    if (s.empty()) return out;

    if (!s.empty() && s.front() == '(' && s.back() == ')') {
        s = s.substr(1, s.size() - 2);
    }

    std::vector<double> vals;
    std::string token;
    std::istringstream iss(s);
    while (std::getline(iss, token, ',')) {
        token = trim(token);
        if (!token.empty()) vals.push_back(std::atof(token.c_str()));
    }
    if (!vals.empty()) out.value = vals[0];
    if (vals.size() > 1) out.stat = vals[1];
    if (vals.size() > 2) out.sys  = vals[2];
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

static Triple invert_triple(const Triple &in) {
    Triple out{0.0, 0.0, 0.0};
    if (in.value == 0.0) {
        return out;
    }

    const double inv = 1.0 / in.value;
    const double factor = std::fabs(1.0 / (in.value * in.value));

    out.value = inv;
    if (in.stat > 0.0) {
        out.stat = factor * in.stat;
    }
    if (in.sys > 0.0) {
        out.sys = factor * in.sys;
    }
    return out;
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
// Luminosity helpers
// -----------------------------------------------------------------------------

static Triple load_rga_lumi_file(const std::string &path,
                                 bool total_from_pos_neg) {
    std::ifstream ifs(path);
    if (!ifs) {
        std::cerr << "[cross_sections] FATAL: cannot open lumi file: "
                  << path << "\n";
        throw std::runtime_error("cannot open lumi file: " + path);
    }

    double sum_total_col = 0.0;
    double sum_pos       = 0.0;
    double sum_neg       = 0.0;

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

        double total = std::atof(fields[1].c_str());
        double pos   = std::atof(fields[2].c_str());
        double neg   = std::atof(fields[3].c_str());

        sum_total_col += total;
        sum_pos       += pos;
        sum_neg       += neg;
        ++n_lines;
    }

    double final_total = total_from_pos_neg ? (sum_pos + sum_neg) : sum_total_col;

    std::cout << "[cross_sections] Loaded lumi from " << path
              << " over " << n_lines << " runs: "
              << "total=" << final_total
              << " pos="   << sum_pos
              << " neg="   << sum_neg << "\n";

    Triple out;
    out.value = final_total;
    out.stat  = sum_pos;
    out.sys   = sum_neg;
    return out;
}

using LumiMap = std::map<std::string, Triple>;

// -----------------------------------------------------------------------------
// Luminosity map construction
// -----------------------------------------------------------------------------

LumiMap build_lumi_map() {
    LumiMap m;

    const std::string base = "imports/integrated_luminosity";

    try {
        m["Fa18 Inb"]      = load_rga_lumi_file(base + "/rga_fa18_inb.txt",  true);
        m["Fa18 Out"]      = load_rga_lumi_file(base + "/rga_fa18_out.txt",  true);
        m["Fa18 Inb Supp"] = Triple{0.0, 0.0, 0.0};
        m["Sp18 Inb"]      = load_rga_lumi_file(base + "/rga_sp18_inb.txt", false);
        m["Sp18 Out"]      = load_rga_lumi_file(base + "/rga_sp18_out.txt", false);
        m["Sp19 Inb"]      = load_rga_lumi_file(base + "/rga_sp19_inb.txt", true);
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

    m["Fa18"]     = sum_labels({"Fa18 Inb", "Fa18 Out", "Fa18 Inb Supp"});
    m["Sp18"]     = sum_labels({"Sp18 Inb", "Sp18 Out"});
    m["10.6 GeV"] = sum_labels({
        "Fa18 Inb", "Fa18 Out", "Fa18 Inb Supp",
        "Sp18 Inb", "Sp18 Out"
    });

    m["10.2 GeV"] = sum_labels({"Sp19 Inb"});

    std::cout << "[cross_sections] build_lumi_map summary:"
              << " Fa18 Inb total=" << m["Fa18 Inb"].value
              << " Sp19 Inb total=" << m["Sp19 Inb"].value
              << " Sp18 Inb total=" << m["Sp18 Inb"].value
              << " (10.6 GeV total=" << m["10.6 GeV"].value
              << ", 10.2 GeV total=" << m["10.2 GeV"].value << ")\n";

    return m;
}

// -----------------------------------------------------------------------------
// Theory model wrappers
// -----------------------------------------------------------------------------

static Helicity helicity_from_string(const std::string &h) {
    if (h == "pos")   return Helicity::Plus;
    if (h == "neg")   return Helicity::Minus;
    return Helicity::Unpol;
}

static double eval_bh_xs(double Ebeam,
                         double xB, double Q2, double t_pos,
                         double phi_rad) {
    const double phi_deg = phi_rad * 180.0 / M_PI;
    return vgg_bh_only(xB, Q2, t_pos, phi_deg, Ebeam);
}

static double eval_km_xs(double Ebeam,
                         double xB, double Q2, double t_pos,
                         double phi_rad, const std::string &h) {
    const double phi_deg = phi_rad * 180.0 / M_PI;
    Helicity hel = helicity_from_string(h);
    return km15_xs(xB, Q2, t_pos, phi_deg, Ebeam, hel);
}

static double eval_vgg_xs(double Ebeam,
                          double xB, double Q2, double t_pos,
                          double phi_rad, const std::string &h) {
    const double phi_deg = phi_rad * 180.0 / M_PI;
    Helicity hel = helicity_from_string(h);
    return vgg_xs(xB, Q2, t_pos, phi_deg, Ebeam, hel);
}

// -----------------------------------------------------------------------------
// JSON helpers
// -----------------------------------------------------------------------------

static void ensure_dir(const fs::path &p) {
    if (!fs::exists(p)) {
        fs::create_directories(p);
    }
}

// -----------------------------------------------------------------------------
// Theory JSON generator
// -----------------------------------------------------------------------------

static bool write_theory_json_for_energy(const std::string &csv_main,
                                         const std::string &theory_json_root,
                                         double Ebeam,
                                         const std::string &energy_label) {
    std::vector<double> phi_deg;
    phi_deg.reserve(14);

    phi_deg.push_back(0.5);
    for (int k = 0; k < 12; ++k) {
        double phi = 15.0 + 30.0 * (double)k;
        phi_deg.push_back(phi);
    }
    phi_deg.push_back(395.5);

    const int n_phi = (int)phi_deg.size();

    std::ifstream ifs(csv_main);
    if (!ifs) {
        std::cerr << "[cross_sections] FATAL: cannot open " << csv_main
                  << " for theory JSON generation (energy " << energy_label
                  << ", Ebeam=" << Ebeam << ").\n";
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
                  << " is empty in write_theory_json_for_energy (energy "
                  << energy_label << ").\n";
        return false;
    }

    std::vector<std::string> header = split_csv_line(lines[0]);

    int c_xb_min = -1, c_xb_max = -1;
    int c_q2_min = -1, c_q2_max = -1;
    int c_t_min  = -1, c_t_max  = -1;
    try {
        c_xb_min = find_col(header, "xBmin");
        c_xb_max = find_col(header, "xBmax");
        c_q2_min = find_col(header, "Q2min");
        c_q2_max = find_col(header, "Q2max");
        c_t_min  = find_col(header, "t_abs_min");
        c_t_max  = find_col(header, "t_abs_max");
    } catch (const std::exception &e) {
        std::cerr << "[cross_sections] FATAL: " << e.what()
                  << " in write_theory_json_for_energy (energy "
                  << energy_label << ").\n";
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
            double frac = 100.0 * (double)row
                          / (double)(n_rows - 1);
            if (frac >= next_pct) {
                std::cout << "[cross_sections] theory JSON ("
                          << energy_label << "): ~" << next_pct
                          << "% of rows processed (row " << row
                          << " / " << (n_rows - 1) << ")\n";
                next_pct += 10;
            }
        }

        std::vector<std::string> fields = split_csv_line(lines[row]);
        if (fields.size() != header.size()) {
            std::cerr << "[cross_sections] WARNING: row " << row
                      << " has " << fields.size() << " fields, expected "
                      << header.size() << " (skipping for theory JSON "
                      << energy_label << ").\n";
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

        std::vector<double> bh_unpol(n_phi), bh_pos(n_phi), bh_neg(n_phi);
        std::vector<double> km_unpol(n_phi), km_pos(n_phi), km_neg(n_phi);
        std::vector<double> vgg_unpol(n_phi), vgg_pos(n_phi), vgg_neg(n_phi);

        for (int i = 0; i < n_phi; ++i) {
            double phi_d = phi_deg[i];
            double phi_r = phi_d * M_PI / 180.0;

            double bh = eval_bh_xs(Ebeam, xB_mid, Q2_mid, t_mid, phi_r);
            bh_unpol[i] = bh;
            bh_pos[i]   = bh;
            bh_neg[i]   = bh;

            km_unpol[i] = eval_km_xs(Ebeam, xB_mid, Q2_mid, t_mid, phi_r, "unpol");
            km_pos[i]   = eval_km_xs(Ebeam, xB_mid, Q2_mid, t_mid, phi_r, "pos");
            km_neg[i]   = eval_km_xs(Ebeam, xB_mid, Q2_mid, t_mid, phi_r, "neg");

            vgg_unpol[i] = eval_vgg_xs(Ebeam, xB_mid, Q2_mid, t_mid, phi_r, "unpol");
            vgg_pos[i]   = eval_vgg_xs(Ebeam, xB_mid, Q2_mid, t_mid, phi_r, "pos");
            vgg_neg[i]   = eval_vgg_xs(Ebeam, xB_mid, Q2_mid, t_mid, phi_r, "neg");
        }

        json bh_json, km_json, vgg_json;
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
                  << " for writing theory JSON (energy "
                  << energy_label << ").\n";
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
        std::cerr << "[cross_sections] ERROR: regenerate_theory_jsons encountered "
                  << "a failure (10.6 ok=" << ok_106
                  << ", 10.2 ok=" << ok_102 << ").\n";
        return false;
    }

    std::cout << "[cross_sections] regenerate_theory_jsons completed successfully.\n";
    return true;
}

// -----------------------------------------------------------------------------
// Cross section computation
// -----------------------------------------------------------------------------

bool compute_cross_sections(const std::string &csv_main,
                            const LumiMap &lumi_map) {
    std::ifstream ifs(csv_main);
    if (!ifs) {
        std::cerr << "[cross_sections] ERROR: cannot open " << csv_main
                  << " for reading.\n";
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
                  << " is empty.\n";
        return false;
    }

    std::vector<std::string> header = split_csv_line(lines[0]);
    const size_t n_data_rows = (lines.size() > 0) ? (lines.size() - 1) : 0;

    LeeFradFbinTable lee_table;
    if (kUseLeeFradFbin) {
        try {
            lee_table = load_lee_frad_fbin_table(kLeeAllBinCsvPath);
        } catch (const std::exception &e) {
            std::cerr << "[cross_sections] FATAL: failed to load Lee Frad/Fbin "
                      << "table from " << kLeeAllBinCsvPath << ": "
                      << e.what() << "\n";
            return false;
        }

        if (lee_table.size() < n_data_rows) {
            std::cerr << "[cross_sections] FATAL: Lee Frad/Fbin table has fewer rows ("
                      << lee_table.size()
                      << ") than dvcs_pass2_analysis.csv (" << n_data_rows
                      << " data rows). Need at least one row per data bin.\n";
            return false;
        }

        if (lee_table.size() > n_data_rows) {
            std::cout << "[cross_sections] INFO: Lee Frad/Fbin table has "
                      << lee_table.size()
                      << " rows while " << csv_main
                      << " has " << n_data_rows
                      << " data rows. Assuming the first "
                      << n_data_rows
                      << " rows correspond one-to-one to valid bins and "
                      << "ignoring the extra trailing rows (valid_bin==0).\n";
        }
    }

    int c_xb_min  = -1, c_xb_max  = -1;
    int c_q2_min  = -1, c_q2_max  = -1;
    int c_t_min   = -1, c_t_max   = -1;
    int c_phi_min = -1, c_phi_max = -1;
    try {
        c_xb_min  = find_col(header, "xBmin");
        c_xb_max  = find_col(header, "xBmax");
        c_q2_min  = find_col(header, "Q2min");
        c_q2_max  = find_col(header, "Q2max");
        c_t_min   = find_col(header, "t_abs_min");
        c_t_max   = find_col(header, "t_abs_max");
        c_phi_min = find_col(header, "phimin");
        c_phi_max = find_col(header, "phimax");
    } catch (const std::exception &e) {
        std::cerr << "[cross_sections] FATAL: " << e.what() << "\n";
        return false;
    }

    int c_xb_idx = find_col_optional(header, "xB index");

    int c_frad_new_106 = find_col_optional(header, "Frad, 10.6 GeV");
    int c_vbin_106     = find_col_optional(header, "bin_volume, 10.6 GeV");
    int c_fbin_new_106 = find_col_optional(header, "Fbin, 10.6 GeV");
    int c_frad_new_102 = find_col_optional(header, "Frad, 10.2 GeV");
    int c_vbin_102     = find_col_optional(header, "bin_volume, 10.2 GeV");
    int c_fbin_new_102 = find_col_optional(header, "Fbin, 10.2 GeV");

    if (c_vbin_106 < 0) {
        std::cerr << "[cross_sections] FATAL: Missing required column "
                  << "\"bin_volume, 10.6 GeV\".\n";
        return false;
    }
    if (c_vbin_102 < 0) {
        std::cerr << "[cross_sections] FATAL: Missing required column "
                  << "\"bin_volume, 10.2 GeV\".\n";
        return false;
    }

    const std::vector<std::string> lumi_col_names_required = {
        "integrated luminosity, Fa18 Inb (nC)",
        "integrated luminosity, Fa18 Out (nC)",
        "integrated luminosity, Sp19 Inb (nC)",
        "integrated luminosity, Sp18 Inb (nC)",
        "integrated luminosity, Sp18 Out (nC)",
        "integrated luminosity, Fa18 (nC)",
        "integrated luminosity, Sp18 (nC)",
        "integrated luminosity, 10.6 GeV (nC)"
    };
    const std::string lumi_col_fa18_supp =
        "integrated luminosity, Fa18 Inb Supp (nC)";
    const std::string lumi_col_10p2 =
        "integrated luminosity, 10.2 GeV (nC)";

    std::map<std::string, int> lumi_col_idx;
    try {
        for (const auto &col_name : lumi_col_names_required) {
            int idx = find_col(header, col_name);
            lumi_col_idx[col_name] = idx;
        }
    } catch (const std::exception &e) {
        std::cerr << "[cross_sections] FATAL: " << e.what() << "\n";
        return false;
    }

    int idx_fa18_supp = find_col_optional(header, lumi_col_fa18_supp);
    if (idx_fa18_supp >= 0) {
        lumi_col_idx[lumi_col_fa18_supp] = idx_fa18_supp;
    }
    int idx_10p2 = find_col_optional(header, lumi_col_10p2);
    if (idx_10p2 >= 0) {
        lumi_col_idx[lumi_col_10p2] = idx_10p2;
    }

    const std::vector<std::string> helicities = {"unpol", "pos", "neg"};

    const std::vector<std::string> labels = {
        "Fa18 Inb", "Fa18 Out", "Fa18 Inb Supp",
        "Sp18 Inb", "Sp18 Out", "Sp19 Inb",
        "Fa18", "Sp18", "10.6 GeV", "10.2 GeV"
    };

    struct ColPair {
        int yield_idx;
        int xs_idx;
    };
    std::map<std::string, ColPair> colmap;

    for (const auto &L : labels) {
        if (L == "10.2 GeV") {
            continue;
        }

        const std::string YL = yield_label_for(L);
        for (const auto &h : helicities) {
            std::string key = L + "|" + h;

            std::string yield_col =
                "acceptance corrected yield, ep->epg, exp, " + YL + ", " + h;
            std::string xs_col =
                "cross sections, ep->epg, exp, " + L + ", " + h;

            int iy = find_col_optional(header, yield_col);
            int ix = find_col_optional(header, xs_col);

            if (iy >= 0 && ix >= 0) {
                colmap[key] = {iy, ix};
            } else if (iy >= 0 && ix < 0) {
                std::cerr << "[cross_sections] FATAL: column \"" << xs_col
                          << "\" is missing while \"" << yield_col
                          << "\" exists. Please add the cross sections column.\n";
                return false;
            }
        }
    }

    bool need_106 = false;
    bool need_102 = false;
    for (const auto &kv : colmap) {
        const std::string &key = kv.first;
        std::string label = key.substr(0, key.find('|'));
        double E = beam_energy_for_label(label);
        if (E > 10.3) need_106 = true;
        else          need_102 = true;
    }

    if (!kUseLeeFradFbin) {
        if (need_106 && c_frad_new_106 < 0) {
            std::cerr << "[cross_sections] FATAL: Missing required column "
                      << "\"Frad, 10.6 GeV\" while kUseLeeFradFbin is false.\n";
            return false;
        }
        if (need_106 && c_fbin_new_106 < 0) {
            std::cerr << "[cross_sections] FATAL: Missing required column "
                      << "\"Fbin, 10.6 GeV\" while kUseLeeFradFbin is false.\n";
            return false;
        }
        if (need_102 && c_frad_new_102 < 0) {
            std::cerr << "[cross_sections] FATAL: Missing required column "
                      << "\"Frad, 10.2 GeV\" while kUseLeeFradFbin is false.\n";
            return false;
        }
        if (need_102 && c_fbin_new_102 < 0) {
            std::cerr << "[cross_sections] FATAL: Missing required column "
                      << "\"Fbin, 10.2 GeV\" while kUseLeeFradFbin is false.\n";
            return false;
        }
    }

    std::vector<std::string> out_lines;
    out_lines.reserve(lines.size());
    out_lines.push_back(lines[0]);

    std::cout << "[cross_sections] compute_cross_sections: total data rows = "
              << n_data_rows << " (excluding header)\n";

    int next_pct = 1;

    for (size_t row = 1; row < lines.size(); ++row) {
        if (lines[row].empty()) {
            out_lines.push_back(lines[row]);
            continue;
        }

        if (n_data_rows > 0 && next_pct <= 100) {
            double frac = 100.0 * (double)row /
                          (double)n_data_rows;
            if (frac >= next_pct) {
                std::cout << "[cross_sections] compute_cross_sections: ~"
                          << next_pct << "% of rows processed (row "
                          << row << " / " << n_data_rows << ")\n";
                next_pct += 10;
            }
        }

        std::vector<std::string> fields = split_csv_line(lines[row]);
        if (fields.size() != header.size()) {
            std::cerr << "[cross_sections] WARNING: row " << row
                      << " has " << fields.size() << " fields, expected "
                      << header.size() << " - copying unchanged.\n";
            out_lines.push_back(lines[row]);
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

        (void)xB_mid;
        (void)Q2_mid;
        (void)t_mid;

        Triple frad_lee_eff{0.0, 0.0, 0.0};
        Triple fbin_lee_eff{1.0, 0.0, 0.0};
        Triple frad_106{0.0, 0.0, 0.0};
        Triple frad_102{0.0, 0.0, 0.0};
        Triple fbin_106{1.0, 0.0, 0.0};
        Triple fbin_102{1.0, 0.0, 0.0};
        Triple vbin_106{0.0, 0.0, 0.0};
        Triple vbin_102{0.0, 0.0, 0.0};

        if (c_vbin_106 >= 0) {
            vbin_106 = parse_tuple3(fields[c_vbin_106]);
        }
        if (c_vbin_102 >= 0) {
            vbin_102 = parse_tuple3(fields[c_vbin_102]);
        }

        if (kUseLeeFradFbin) {
            const LeeFradFbinRow &lr = lee_table[row - 1];

            if (lr.frad == 0.0) {
                std::cerr << "[cross_sections] FATAL: Lee Frad is zero at row "
                          << row << " in " << kLeeAllBinCsvPath
                          << ". This makes Frad ill-defined.\n";
                return false;
            }
            if (lr.fbin == 0.0) {
                std::cerr << "[cross_sections] FATAL: Lee Fbin is zero at row "
                          << row << " in " << kLeeAllBinCsvPath
                          << ". This makes Fbin ill-defined.\n";
                return false;
            }

            frad_lee_eff.value = 1.0 / lr.frad;
            frad_lee_eff.stat  = 0.0;
            frad_lee_eff.sys   = 0.0;

            fbin_lee_eff.value = 1.0 / lr.fbin;
            fbin_lee_eff.stat  = 0.0;
            fbin_lee_eff.sys   = 0.0;
        } else {
            if (c_frad_new_106 >= 0) {
                frad_106 = parse_tuple3(fields[c_frad_new_106]);
            }
            if (c_frad_new_102 >= 0) {
                frad_102 = parse_tuple3(fields[c_frad_new_102]);
            }
            if (c_fbin_new_106 >= 0) {
                fbin_106 = parse_tuple3(fields[c_fbin_new_106]);
            }
            if (c_fbin_new_102 >= 0) {
                fbin_102 = parse_tuple3(fields[c_fbin_new_102]);
            }
        }

        for (const auto &L : labels) {
            auto lumi_it = lumi_map.find(L);
            if (lumi_it == lumi_map.end()) continue;
            const Triple &Lumi = lumi_it->second;

            std::string lumi_col_name;
            if (L == "Fa18 Inb")      lumi_col_name = "integrated luminosity, Fa18 Inb (nC)";
            else if (L == "Fa18 Out") lumi_col_name = "integrated luminosity, Fa18 Out (nC)";
            else if (L == "Fa18 Inb Supp") lumi_col_name = "integrated luminosity, Fa18 Inb Supp (nC)";
            else if (L == "Sp19 Inb") lumi_col_name = "integrated luminosity, Sp19 Inb (nC)";
            else if (L == "Sp18 Inb") lumi_col_name = "integrated luminosity, Sp18 Inb (nC)";
            else if (L == "Sp18 Out") lumi_col_name = "integrated luminosity, Sp18 Out (nC)";
            else if (L == "Fa18")     lumi_col_name = "integrated luminosity, Fa18 (nC)";
            else if (L == "Sp18")     lumi_col_name = "integrated luminosity, Sp18 (nC)";
            else if (L == "10.6 GeV") lumi_col_name = "integrated luminosity, 10.6 GeV (nC)";
            else if (L == "10.2 GeV") lumi_col_name = "integrated luminosity, 10.2 GeV (nC)";

            auto it_lumi_col = lumi_col_idx.find(lumi_col_name);
            if (it_lumi_col != lumi_col_idx.end()) {
                fields[it_lumi_col->second] =
                    tuple3_to_cell(Lumi.value, Lumi.stat, Lumi.sys);
            }

            for (const auto &h : helicities) {
                std::string key = L + "|" + h;
                auto it = colmap.find(key);
                if (it == colmap.end()) continue;

                int iy = it->second.yield_idx;
                int ix = it->second.xs_idx;
                if (iy < 0 || ix < 0) continue;

                Triple Y = parse_tuple3(fields[iy]);
                if (Y.value <= 0.0) continue;

                double lumi_val = 0.0;
                if (h == "unpol")      lumi_val = Lumi.value;
                else if (h == "pos")   lumi_val = Lumi.stat;
                else if (h == "neg")   lumi_val = Lumi.sys;

                if (lumi_val <= 0.0) {
                    std::cerr << "[cross_sections] WARNING: zero or negative "
                              << "lumi for label=" << L
                              << " helicity=" << h
                              << " at row " << row << ".\n";
                    continue;
                }

                double E = beam_energy_for_label(L);
                const Triple *vbin_ptr  = nullptr;
                const Triple *frad_ptr  = nullptr;
                Triple        fbin_eff{1.0, 0.0, 0.0};

                if (E > 10.3) {
                    vbin_ptr = &vbin_106;
                    if (kUseLeeFradFbin) {
                        frad_ptr = &frad_lee_eff;
                        fbin_eff = fbin_lee_eff;
                    } else {
                        frad_ptr = &frad_106;
                        fbin_eff = fbin_106;
                    }
                } else {
                    vbin_ptr = &vbin_102;
                    if (kUseLeeFradFbin) {
                        frad_ptr = &frad_lee_eff;
                        fbin_eff = fbin_lee_eff;
                    } else {
                        frad_ptr = &frad_102;
                        fbin_eff = fbin_102;
                    }
                }

                if (!vbin_ptr || !frad_ptr) {
                    std::cerr << "[cross_sections] WARNING: missing vbin/frad "
                              << "for label=" << L
                              << " helicity=" << h
                              << " at row " << row << ".\n";
                    continue;
                }

                if (vbin_ptr->value <= 0.0) {
                    std::cerr << "[cross_sections] WARNING: zero or negative "
                              << "bin volume at row " << row
                              << " for label=" << L
                              << " helicity=" << h << ".\n";
                    continue;
                }
                if (frad_ptr->value <= 0.0) {
                    std::cerr << "[cross_sections] WARNING: zero or negative "
                              << "Frad at row " << row
                              << " for label=" << L
                              << " helicity=" << h << ".\n";
                    continue;
                }

                double denom = vbin_ptr->value * lumi_val;
                if (denom <= 0.0) {
                    std::cerr << "[cross_sections] WARNING: zero or negative denom "
                              << "(vbin * lumi) at row " << row
                              << " for label=" << L << ", helicity=" << h << ".\n";
                    continue;
                }

                double sigma = Y.value * frad_ptr->value * fbin_eff.value / denom;

                double rel2 = 0.0;
                if (Y.value > 0.0 && Y.stat > 0.0) {
                    double r = Y.stat / Y.value;
                    rel2 += r * r;
                }
                if (frad_ptr->value > 0.0 && frad_ptr->stat > 0.0) {
                    double r = frad_ptr->stat / frad_ptr->value;
                    rel2 += r * r;
                }
                if (fbin_eff.value > 0.0 && fbin_eff.stat > 0.0) {
                    double r = fbin_eff.stat / fbin_eff.value;
                    rel2 += r * r;
                }

                double sigma_stat = (rel2 > 0.0) ? sigma * std::sqrt(rel2) : 0.0;
                double sigma_sys  = 0.0;

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
        std::cerr << "[cross_sections] ERROR: cannot open " << tmp_path
                  << " for writing.\n";
        return false;
    }
    for (const auto &lout : out_lines) {
        ofs << lout << "\n";
    }
    ofs.close();

    fs::rename(tmp_path, csv_path);
    std::cout << "[cross_sections] Updated CSV with cross sections and luminosities: "
              << csv_main << "\n";

    return true;
}

// -----------------------------------------------------------------------------
// Plotting
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
    bool   have_theory_row = false;
};

using QTKey = std::pair<Range,Range>;

struct XSGroupByXB {
    std::map<QTKey, BinData> bins;
    int xb_index = -1;
};

struct TheoryCurves {
    std::vector<double> phi_deg;
    std::vector<double> bh_unpol, bh_pos, bh_neg;
    std::vector<double> km_unpol, km_pos, km_neg;
    std::vector<double> vgg_unpol, vgg_pos, vgg_neg;
};

static std::string theory_energy_label_for(const std::string &label) {
    if (label == "Sp19 Inb" || label == "10.2 GeV") {
        return "10.2 GeV";
    }
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
        std::cerr << "[cross_sections] WARNING: no theory JSON file for label \""
                  << label << "\" (energy \"" << energy_label
                  << "\") at " << file.string() << "\n";
        return out;
    }

    std::ifstream ifs(file);
    if (!ifs) {
        std::cerr << "[cross_sections] WARNING: failed to open theory JSON for label \""
                  << label << "\" (energy \"" << energy_label
                  << "\") at " << file.string() << "\n";
        return out;
    }

    json j;
    try {
        ifs >> j;
    } catch (...) {
        std::cerr << "[cross_sections] WARNING: malformed theory JSON for label \""
                  << label << "\" (energy \"" << energy_label
                  << "\") at " << file.string() << "\n";
        return out;
    }

    std::vector<double> phi_deg = j.value("phi_deg", std::vector<double>{});
    if (phi_deg.empty()) {
        std::cerr << "[cross_sections] WARNING: theory JSON for label \""
                  << label << "\" (energy \"" << energy_label
                  << "\") has empty phi_deg.\n";
        return out;
    }

    if (!j.contains("rows") || !j["rows"].is_object()) {
        std::cerr << "[cross_sections] WARNING: theory JSON for label \""
                  << label << "\" (energy \"" << energy_label
                  << "\") has no 'rows' object.\n";
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
              << "\" using energy \"" << energy_label << "\" from "
              << file.string()
              << " with " << out.size() << " rows.\n";

    return out;
}

// -----------------------------------------------------------------------------
// View mode for canvases
// -----------------------------------------------------------------------------

enum class XSecPanelMode {
    All,
    UnpolOnly,
    PosOnly,
    NegOnly
};

static std::pair<double,double> compute_yrange_for_bin(
    const BinData *bin,
    const std::map<size_t, TheoryCurves> &theory,
    XSecPanelMode mode)
{
    double ymin = std::numeric_limits<double>::max();
    double ymax = 0.0;

    auto update_from_points = [&](const std::vector<Point> &v) {
        for (const auto &p : v) {
            if (p.xs > 0.0) {
                double y_low  = std::max(1e-12, p.xs - p.xs_err);
                double y_high = p.xs + p.xs_err;
                if (y_low  > 0.0 && y_low  < ymin) ymin = y_low;
                if (y_high > ymax) ymax = y_high;
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
        double log_min = std::pow(10.0, std::floor(std::log10(ymin)));
        double log_max = std::pow(10.0, std::ceil(std::log10(ymax)));
        ymin = std::max(1e-4, log_min);
        ymax = log_max;
    }

    return std::make_pair(ymin, 1.2 * ymax);
}

static TGraphErrors* make_xsec_graph(const std::vector<Point> &v,
                                     int mstyle, int mcolor) {
    if (v.empty()) return nullptr;
    int N = (int)v.size();
    std::vector<double> x(N), y(N), ex(N), ey(N);
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
// Canvas builder
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
    int nPads)
{
    const auto &bins_for_xB = group.bins;
    if (bins_for_xB.empty()) return;

    int W = 280 * ncols + 160;
    int H = 260 * nrows + 260;

    const int MIN_W = 1200;
    const int MIN_H = 900;
    if (W < MIN_W) W = MIN_W;
    if (H < MIN_H) H = MIN_H;

    double titleSize = 0.18;
    double legendTextSize = 0.11;
    double cellLabelSize = 0.070;

    if (nPads <= 4) {
        titleSize = 0.14;
        legendTextSize = 0.09;
        cellLabelSize = 0.060;
    }
    if (nPads == 1) {
        titleSize = 0.12;
        legendTextSize = 0.085;
        cellLabelSize = 0.055;
    }

    titleSize *= 0.5;
    legendTextSize *= 0.5;

    std::ostringstream cname;
    cname << "c_xsec_";
    if (mode == XSecPanelMode::UnpolOnly) cname << "unpol_";
    else if (mode == XSecPanelMode::PosOnly) cname << "pos_";
    else if (mode == XSecPanelMode::NegOnly) cname << "neg_";
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

    TGraphErrors dummy_unpol, dummy_pos, dummy_neg;
    dummy_unpol.SetMarkerStyle(20);
    dummy_unpol.SetLineWidth(2);
    dummy_unpol.SetMarkerColor(kBlack);
    dummy_unpol.SetLineColor(kBlack);

    dummy_pos.SetMarkerStyle(24);
    dummy_pos.SetLineWidth(2);
    dummy_pos.SetMarkerColor(kRed+1);
    dummy_pos.SetLineColor(kRed+1);

    dummy_neg.SetMarkerStyle(25);
    dummy_neg.SetLineWidth(2);
    dummy_neg.SetMarkerColor(kBlue+1);
    dummy_neg.SetLineColor(kBlue+1);

    TGraph dummy_bh, dummy_km_unpol, dummy_km_pos, dummy_km_neg;
    TGraph dummy_vgg_unpol, dummy_vgg_pos, dummy_vgg_neg;

    dummy_bh.SetLineWidth(2);
    dummy_bh.SetLineStyle(2);
    dummy_bh.SetLineColor(kGreen+2);

    dummy_km_unpol.SetLineWidth(2);
    dummy_km_unpol.SetLineStyle(1);
    dummy_km_unpol.SetLineColor(kMagenta+1);

    dummy_km_pos.SetLineWidth(2);
    dummy_km_pos.SetLineStyle(2);
    dummy_km_pos.SetLineColor(kMagenta+1);

    dummy_km_neg.SetLineWidth(2);
    dummy_km_neg.SetLineStyle(3);
    dummy_km_neg.SetLineColor(kMagenta+1);

    dummy_vgg_unpol.SetLineWidth(2);
    dummy_vgg_unpol.SetLineStyle(1);
    dummy_vgg_unpol.SetLineColor(kOrange+7);

    dummy_vgg_pos.SetLineWidth(2);
    dummy_vgg_pos.SetLineStyle(2);
    dummy_vgg_pos.SetLineColor(kOrange+7);

    dummy_vgg_neg.SetLineWidth(2);
    dummy_vgg_neg.SetLineStyle(3);
    dummy_vgg_neg.SetLineColor(kOrange+7);

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

        legKM->AddEntry(&dummy_bh,       "BH unpolarized",  "l");
        legKM->AddEntry(&dummy_km_unpol, "KM unpolarized",  "l");
        legKM->AddEntry(&dummy_km_pos,   "KM + helicity",   "l");
        legKM->AddEntry(&dummy_km_neg,   "KM - helicity",   "l");

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

        legKM->AddEntry(&dummy_bh,     "BH",             "l");
        legKM->AddEntry(&dummy_km_pos, "KM + helicity",  "l");

        legVGG->AddEntry(&dummy_vgg_pos, "VGG + helicity", "l");
    } else if (mode == XSecPanelMode::NegOnly) {
        legData->AddEntry(&dummy_neg, "data - helicity", "lep");

        legKM->AddEntry(&dummy_bh,     "BH",             "l");
        legKM->AddEntry(&dummy_km_neg, "KM - helicity",  "l");

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

            if (!bin_ptr) {
                continue;
            }

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
                if (g_unpol) g_unpol->Draw("P SAME");
            }
            if (mode == XSecPanelMode::All || mode == XSecPanelMode::PosOnly) {
                TGraphErrors *g_pos = make_xsec_graph(bin.pos, 24, kRed+1);
                if (g_pos) g_pos->Draw("P SAME");
            }
            if (mode == XSecPanelMode::All || mode == XSecPanelMode::NegOnly) {
                TGraphErrors *g_neg = make_xsec_graph(bin.neg, 25, kBlue+1);
                if (g_neg) g_neg->Draw("P SAME");
            }

            if (bin.have_theory_row) {
                auto it_th = theory.find(bin.theory_row);
                if (it_th != theory.end()) {
                    const TheoryCurves &tc = it_th->second;

                    auto draw_curve = [&](const std::vector<double> &ys,
                                          int lstyle, int lcolor) {
                        if (tc.phi_deg.size() != ys.size() || ys.empty()) return;
                        int M = (int)ys.size();
                        std::vector<double> xp(M), yp(M);
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
                        draw_curve(tc.bh_unpol, 2, kGreen+2);
                        draw_curve(tc.km_unpol, 1, kMagenta+1);
                        draw_curve(tc.km_pos,   2, kMagenta+1);
                        draw_curve(tc.km_neg,   3, kMagenta+1);
                        draw_curve(tc.vgg_unpol, 1, kOrange+7);
                        draw_curve(tc.vgg_pos,   2, kOrange+7);
                        draw_curve(tc.vgg_neg,   3, kOrange+7);
                    } else if (mode == XSecPanelMode::UnpolOnly) {
                        draw_curve(tc.bh_unpol, 2, kGreen+2);
                        draw_curve(tc.km_unpol, 1, kMagenta+1);
                        draw_curve(tc.vgg_unpol, 1, kOrange+7);
                    } else if (mode == XSecPanelMode::PosOnly) {
                        draw_curve(tc.bh_unpol, 2, kGreen+2);
                        draw_curve(tc.km_pos,   2, kMagenta+1);
                        draw_curve(tc.vgg_pos,  2, kOrange+7);
                    } else if (mode == XSecPanelMode::NegOnly) {
                        draw_curve(tc.bh_unpol, 2, kGreen+2);
                        draw_curve(tc.km_neg,   3, kMagenta+1);
                        draw_curve(tc.vgg_neg,  3, kOrange+7);
                    }
                }
            }
        }
    }

    std::ostringstream fname;
    fname << "cross_sections_";
    if (mode == XSecPanelMode::UnpolOnly) fname << "unpol_";
    else if (mode == XSecPanelMode::PosOnly) fname << "pos_";
    else if (mode == XSecPanelMode::NegOnly) fname << "neg_";
    fname << canonical_period_dir(label)
          << "_xB_" << xb_idx_for_name << ".png";

    fs::path outpath = outdir / fname.str();
    c->SaveAs(outpath.string().c_str());

    delete c;
}

// -----------------------------------------------------------------------------
// Main plotting entry point
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

    int c_xb_min = -1, c_xb_max = -1;
    int c_q2_min = -1, c_q2_max = -1;
    int c_t_min  = -1, c_t_max  = -1;
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
    int c_phimin = -1, c_phimax = -1;
    if (c_phiavg < 0) {
        c_phimin = find_col_optional(header, "phimin");
        c_phimax = find_col_optional(header, "phimax");
        if (c_phimin < 0 || c_phimax < 0) {
            std::cerr << "[cross_sections] FATAL: neither \"phiavg, " << label
                      << "\" nor phimin/phimax found in CSV.\n";
            return false;
        }
    }

    int c_xs_unpol = find_col_optional(
        header, "cross sections, ep->epg, exp, " + label + ", unpol");
    int c_xs_pos   = find_col_optional(
        header, "cross sections, ep->epg, exp, " + label + ", pos");
    int c_xs_neg   = find_col_optional(
        header, "cross sections, ep->epg, exp, " + label + ", neg");

    if (c_xs_unpol < 0 && c_xs_pos < 0 && c_xs_neg < 0) {
        std::cerr << "[cross_sections] INFO: no cross section columns found for label "
                  << label << " in CSV; nothing to plot.\n";
        return true;
    }

    if (c_xs_unpol < 0 || c_xs_pos < 0 || c_xs_neg < 0) {
        std::cerr << "[cross_sections] FATAL: cross section columns for label "
                  << label << " are incomplete; expected unpol/pos/neg.\n";
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
        Triple xs_pos   = parse_tuple3(fields[c_xs_pos]);
        Triple xs_neg   = parse_tuple3(fields[c_xs_neg]);

        if (xs_unpol.value <= 0.0 && xs_pos.value <= 0.0 && xs_neg.value <= 0.0) {
            continue;
        }

        XSGroupByXB &group = by_xb[xb_range];

        if (group.xb_index < 0 && c_xb_idx >= 0) {
            group.xb_index = std::atoi(trim(unquote(fields[c_xb_idx])).c_str());
        }

        QTKey key(q2_range, t_range);
        BinData &bin = group.bins[key];

        if (!bin.have_theory_row) {
            bin.theory_row = row;
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
        add_point(xs_pos,   bin.pos);
        add_point(xs_neg,   bin.neg);
    }

    if (by_xb.empty()) {
        std::cerr << "[cross_sections] WARNING: no cross section data found for label "
                  << label << " in CSV (nothing to plot).\n";
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

        std::set<Range> q2_set, t_set;
        for (const auto &kv : bins_for_xB) {
            const QTKey &qt = kv.first;
            q2_set.insert(qt.first);
            t_set.insert(qt.second);
        }
        std::vector<Range> q2_slice(q2_set.begin(), q2_set.end());
        std::vector<Range> t_slice(t_set.begin(), t_set.end());
        if (q2_slice.empty() || t_slice.empty()) continue;

        const int ncols = (int)q2_slice.size();
        const int nrows = (int)t_slice.size();
        const int nPads = ncols * nrows;

        int xb_idx_for_name = (group.xb_index >= 0 ? group.xb_index : xb_canvas_counter);

        make_xsec_canvas_for_mode(label, xb_range, group,
                                  q2_slice, t_slice,
                                  theory, outdir,
                                  xb_idx_for_name,
                                  XSecPanelMode::All,
                                  ncols, nrows, nPads);

        make_xsec_canvas_for_mode(label, xb_range, group,
                                  q2_slice, t_slice,
                                  theory, outdir,
                                  xb_idx_for_name,
                                  XSecPanelMode::UnpolOnly,
                                  ncols, nrows, nPads);

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

        ++xb_canvas_counter;
    }

    std::cout << "[cross_sections] Plotted cross sections (all modes) for label "
              << label << " into " << outdir.string() << "\n";

    return true;
}