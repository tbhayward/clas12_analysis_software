// cross_sections.cpp
// -----------------------------------------------------------------------------
// Cross section calculator and plotting with BH / KM / VGG theory curves.
//
// This module:
//   1) Updates dvcs_pass2_analysis.csv:
//        - Fills integrated luminosity columns for all rows.
//        - Computes cross sections for unpol/pos/neg per label where columns exist.
//        - Stores them as "(value, stat, sys)" with sys = 0.
//        - Propagates stat uncertainties from the acceptance corrected yield
//          and Frad.
//   2) Computes theory curves at bin-mean kinematics for each row for
//      10.6 GeV and 10.2 GeV only, using 12 equally spaced phi points
//      from 0 to 360 degrees, and caches them as JSON under
//      output/jsons/cross_sections/<PeriodDir>/xs_phi_all.json.
//      For each CSV row it stores:
//        - BH unpolarized.
//        - KM: unpol, +, -.
//        - VGG: unpol, +, -.
//   3) For each label, plots xB canvases where:
//        - One canvas per xB bin.
//        - Each canvas is laid out as N_t rows (|t| bins) by N_Q2 columns (Q^2 bins).
//        - Each subpad (a fixed Q^2 and |t| bin) shows:
//            * unpolarized, + helicity, and - helicity cross sections vs phi
//              (with errors) from the CSV.
//            * BH (unpol), KM (unpol/+/−), VGG (unpol/+/−) curves vs phi
//              from the theory JSON for that label (if present).
//        - Y axis is drawn in log scale.
//        - A single legend is drawn once at the top of the canvas.
//
// This file does not add new columns to the CSV; it only updates existing ones.
// If required columns are missing, it fails with a clear error.
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
// Period directory helper (consistent with make_dirs.cpp conventions)
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

    // Fallback: replace spaces with underscores (should not normally be used)
    std::string out = label;
    std::replace(out.begin(), out.end(), ' ', '_');
    return out;
}

// The acceptance corrected yield columns use exactly these period labels,
// except for the combined 10.6 group which uses "2018 (10.6 GeV)".
static std::string yield_label_for(const std::string &label) {
    if (label == "10.6 GeV") {
        return "2018 (10.6 GeV)";
    }
    if (label == "10.2 GeV") {
        // For now this is just Sp19 Inb; if you later add a separate combined label
        // you can adjust this mapping.
        return "Sp19 Inb";
    }
    return label;
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
    while (b < s.size() && std::isspace(static_cast<unsigned char>(s[b]))) ++b;
    size_t e = s.size();
    while (e > b && std::isspace(static_cast<unsigned char>(s[e - 1]))) --e;
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
        if (c == ',' || c == '"' || std::isspace(static_cast<unsigned char>(c))) {
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

// Parse "(value, stat, sys)" -> Triple.
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

// Find column index; fatal if missing.
static int find_col(const std::vector<std::string> &header,
                    const std::string &target) {
    for (size_t i = 0; i < header.size(); ++i) {
        if (trim(unquote(header[i])) == target) {
            return static_cast<int>(i);
        }
    }
    throw std::runtime_error("Missing required column: \"" + target + "\"");
}

// Try to find column; return -1 if missing.
static int find_col_optional(const std::vector<std::string> &header,
                             const std::string &target) {
    for (size_t i = 0; i < header.size(); ++i) {
        if (trim(unquote(header[i])) == target) {
            return static_cast<int>(i);
        }
    }
    return -1;
}

// -----------------------------------------------------------------------------
// Luminosity helpers
// -----------------------------------------------------------------------------

// Helper: load one RGA integrated-luminosity file
//
// Expected line format (comma-separated, no header, one line per run):
//   run, total, pos, neg, ...
//
// We:
//   - ignore empty lines and lines starting with '#'
//   - sum column 2 (total), column 3 (pos), column 4 (neg)
//   - store them as Triple{total, pos, neg}
static Triple load_rga_lumi_file(const std::string &path) {
    std::ifstream ifs(path);
    if (!ifs) {
        std::cerr << "[cross_sections] FATAL: cannot open lumi file: "
                  << path << "\n";
        throw std::runtime_error("cannot open lumi file: " + path);
    }

    double sum_total = 0.0;
    double sum_pos   = 0.0;
    double sum_neg   = 0.0;

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

        double total = std::atof(fields[1].c_str()); // column 2
        double pos   = std::atof(fields[2].c_str()); // column 3
        double neg   = std::atof(fields[3].c_str()); // column 4

        sum_total += total;
        sum_pos   += pos;
        sum_neg   += neg;
        ++n_lines;
    }

    std::cout << "[cross_sections] Loaded lumi from " << path
              << " over " << n_lines << " runs: "
              << "total=" << sum_total
              << " pos="   << sum_pos
              << " neg="   << sum_neg << "\n";

    Triple out;
    out.value = sum_total; // total unpolarized charge
    out.stat  = sum_pos;   // + helicity charge
    out.sys   = sum_neg;   // - helicity charge
    return out;
}

// -----------------------------------------------------------------------------
// Luminosity map construction
//
// For each label L, we store:
//   Lumi.value = total accumulated charge (unpolarized, column 2)
//   Lumi.stat  = + helicity accumulated charge (column 3)
//   Lumi.sys   = - helicity accumulated charge (column 4)
//
// The RGA text files live in imports/integrated_luminosity/:
//
//   rga_fa18_inb.txt   -> "Fa18 Inb"
//   rga_fa18_out.txt   -> "Fa18 Out"
//   rga_sp18_out.txt   -> "Sp18 Out"
//   rga_sp19_inb.txt   -> "Sp19 Inb" (10.2 GeV period)
//
// Combined entries:
//   "Fa18"     = Fa18 Inb + Fa18 Out + Fa18 Inb Supp
//   "Sp18"     = Sp18 Inb + Sp18 Out
//   "10.6 GeV" = all 2018 10.6 periods (Fa18 + Sp18)
//   "10.2 GeV" = Sp19 Inb only
// -----------------------------------------------------------------------------

LumiMap build_lumi_map() {
    LumiMap m;

    const std::string base = "imports/integrated_luminosity";

    try {
        // Fall 2018
        m["Fa18 Inb"]      = load_rga_lumi_file(base + "/rga_fa18_inb.txt");
        m["Fa18 Out"]      = load_rga_lumi_file(base + "/rga_fa18_out.txt");
        // If you have a separate file for Fa18 Inb Supp, load it here.
        // For now we set it to zero and keep it explicit.
        m["Fa18 Inb Supp"] = Triple{0.0, 0.0, 0.0};

        // Spring 2018
        // If you have rga_sp18_inb.txt, load it here. For now set to zero.
        m["Sp18 Inb"]      = Triple{0.0, 0.0, 0.0};
        m["Sp18 Out"]      = load_rga_lumi_file(base + "/rga_sp18_out.txt");

        // Spring 2019 (10.2 GeV)
        m["Sp19 Inb"]      = load_rga_lumi_file(base + "/rga_sp19_inb.txt");
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

    // Combined 2018 periods (10.6 GeV)
    m["Fa18"]     = sum_labels({"Fa18 Inb", "Fa18 Out", "Fa18 Inb Supp"});
    m["Sp18"]     = sum_labels({"Sp18 Inb", "Sp18 Out"});
    m["10.6 GeV"] = sum_labels({
        "Fa18 Inb", "Fa18 Out", "Fa18 Inb Supp",
        "Sp18 Inb", "Sp18 Out"
    });

    // 10.2 GeV = Sp19 Inb only
    m["10.2 GeV"] = sum_labels({"Sp19 Inb"});

    std::cout << "[cross_sections] build_lumi_map summary:"
              << " Fa18 Inb total=" << m["Fa18 Inb"].value
              << " Sp19 Inb total=" << m["Sp19 Inb"].value
              << " (10.6 GeV total=" << m["10.6 GeV"].value
              << ", 10.2 GeV total=" << m["10.2 GeV"].value << ")\n";

    return m;
}

// -----------------------------------------------------------------------------
// Theory model wrappers (BH / KM / VGG)
// -----------------------------------------------------------------------------

static Helicity helicity_from_string(const std::string &h) {
    if (h == "pos")   return Helicity::Plus;
    if (h == "neg")   return Helicity::Minus;
    return Helicity::Unpol;
}

// All wrappers take phi in radians and convert internally to degrees.

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
// Cross section computation (CSV update + theory JSON write)
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

    // Required kinematic columns (for computing bin means)
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

    // Radiative and bin-volume columns (10.6 GeV currently)
    int c_frad = -1;
    int c_vbin = -1;
    try {
        c_frad = find_col(header, "Frad, 10.6 GeV");
        c_vbin = find_col(header, "bin_volume, 10.6 GeV");
    } catch (const std::exception &e) {
        std::cerr << "[cross_sections] FATAL: " << e.what() << "\n";
        return false;
    }

    // Integrated luminosity columns
    // These must already exist in the CSV; we only fill them.
    const std::vector<std::string> lumi_col_names_required = {
        "integrated luminosity, Fa18 Inb (nC)",
        "integrated luminosity, Fa18 Out (nC)",
        "integrated luminosity, Sp19 Inb (nC)",
        "integrated luminosity, Sp18 Inb (nC)",
        "integrated luminosity, Sp18 Out (nC)",
        "integrated luminosity, Fa18 (nC)",
        "integrated luminosity, Sp18 (nC)",
        "integrated luminosity, 10.6 GeV (nC)",
        "integrated luminosity, 10.2 GeV (nC)"
    };
    const std::string lumi_col_fa18_supp =
        "integrated luminosity, Fa18 Inb Supp (nC)";

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
    // Optional Fa18 Inb Supp column
    int idx_fa18_supp = find_col_optional(header, lumi_col_fa18_supp);
    if (idx_fa18_supp >= 0) {
        lumi_col_idx[lumi_col_fa18_supp] = idx_fa18_supp;
    }

    // Helicity labels used in column names
    const std::vector<std::string> helicities = {"unpol", "pos", "neg"};

    // Labels (periods and combined groups) we care about
    const std::vector<std::string> labels = {
        "Fa18 Inb", "Fa18 Out", "Fa18 Inb Supp",
        "Sp18 Inb", "Sp18 Out", "Sp19 Inb",
        "Fa18", "Sp18", "10.6 GeV", "10.2 GeV"
    };

    // Map of columns for each label+helicity:
    //   acceptance corrected yield, ep->epg, exp, <yield_label>, <helicity>
    //   cross sections, ep->epg, exp, <label>, <helicity>
    struct ColPair {
        int yield_idx;
        int xs_idx;
    };
    std::map<std::string, ColPair> colmap;  // key = label + "|" + helicity

    for (const auto &L : labels) {
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
            // If neither exists, we simply do not compute that combination.
        }
    }

    // Theory JSON root
    const std::string theory_root = "output/jsons/cross_sections";
    ensure_dir(theory_root);

    // Labels for which we actually evaluate BH/KM/VGG
    const std::vector<std::string> theory_labels = {
        "10.6 GeV", "10.2 GeV"
    };

    // Precompute phi grid (12 points from 0 to 360 deg)
    const int Nphi = 12;
    std::vector<double> phi_deg(Nphi);
    for (int i = 0; i < Nphi; ++i) {
        double frac = (i + 0.5) / static_cast<double>(Nphi);
        phi_deg[i] = 360.0 * frac;
    }

    // Per-label theory JSON accumulator
    std::map<std::string, json> theory_json;
    for (const auto &L : theory_labels) {
        double Ebeam = (L == "10.2 GeV") ? 10.2 : 10.6;
        json j;
        j["label"]   = L;
        j["Ebeam"]   = Ebeam;
        j["phi_deg"] = phi_deg;
        j["rows"]    = json::object();
        theory_json[L] = j;
    }

    // Output lines
    std::vector<std::string> out_lines;
    out_lines.reserve(lines.size());
    out_lines.push_back(lines[0]);  // header

    const size_t n_data_rows = (lines.size() > 0) ? (lines.size() - 1) : 0;
    std::cout << "[cross_sections] compute_cross_sections: total data rows = "
              << n_data_rows << " (excluding header)\n";

    int next_pct = 10;

    // Process rows
    for (size_t row = 1; row < lines.size(); ++row) {
        if (lines[row].empty()) {
            out_lines.push_back(lines[row]);
            continue;
        }

        // Progress print ~every 10%
        if (n_data_rows > 0 && next_pct <= 100) {
            double frac = 100.0 * static_cast<double>(row) /
                          static_cast<double>(n_data_rows);
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

        // Bin means for theory
        double xbmin = std::atof(trim(unquote(fields[c_xb_min])).c_str());
        double xbmax = std::atof(trim(unquote(fields[c_xb_max])).c_str());
        double q2min = std::atof(trim(unquote(fields[c_q2_min])).c_str());
        double q2max = std::atof(trim(unquote(fields[c_q2_max])).c_str());
        double tmin  = std::atof(trim(unquote(fields[c_t_min])).c_str());
        double tmax  = std::atof(trim(unquote(fields[c_t_max])).c_str());

        double xB_mid = 0.5 * (xbmin + xbmax);
        double Q2_mid = 0.5 * (q2min + q2max);
        double t_mid  = 0.5 * (tmin  + tmax);

        // xB index (if present) is used only for naming plots later
        int xb_idx_val = -1;
        if (c_xb_idx >= 0) {
            xb_idx_val = std::atoi(trim(unquote(fields[c_xb_idx])).c_str());
        }
        (void)xb_idx_val;

        // Global factors for this bin
        Triple frad = parse_tuple3(fields[c_frad]);
        Triple vbin = parse_tuple3(fields[c_vbin]);

        // For each label, compute cross sections and fill lumi
        for (const auto &L : labels) {
            auto lumi_it = lumi_map.find(L);
            if (lumi_it == lumi_map.end()) continue;
            const Triple &Lumi = lumi_it->second;

            // Fill integrated luminosity columns for this label if the CSV has them
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

            // Now compute cross sections for each helicity where we have columns.
            // IMPORTANT: use helicity-specific luminosity:
            //   unpol -> Lumi.value (total charge)
            //   pos   -> Lumi.stat  (+ helicity)
            //   neg   -> Lumi.sys   (- helicity)
            for (const auto &h : helicities) {
                std::string key = L + "|" + h;
                auto it = colmap.find(key);
                if (it == colmap.end()) continue;

                int iy = it->second.yield_idx;
                int ix = it->second.xs_idx;
                if (iy < 0 || ix < 0) continue;

                Triple Y = parse_tuple3(fields[iy]);
                if (Y.value <= 0.0) {
                    // Leave cross section cell as-is (likely empty)
                    continue;
                }

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

                double denom = vbin.value * lumi_val;
                if (denom <= 0.0) {
                    std::cerr << "[cross_sections] WARNING: zero or negative denom "
                              << "(vbin * lumi) at row " << row
                              << " for label=" << L << ", helicity=" << h << ".\n";
                    continue;
                }

                // sigma = yield * Frad / (bin_volume * Lumi_h)
                double sigma = Y.value * frad.value / denom;

                // Stat error from Y and Frad
                double rel2 = 0.0;
                if (Y.value > 0.0 && Y.stat > 0.0) {
                    double r = Y.stat / Y.value;
                    rel2 += r * r;
                }
                if (frad.value > 0.0 && frad.stat > 0.0) {
                    double r = frad.stat / frad.value;
                    rel2 += r * r;
                }
                double sigma_stat = (rel2 > 0.0) ? sigma * std::sqrt(rel2) : 0.0;

                // Sys set to 0 for now; you will edit CSV later.
                double sigma_sys = 0.0;

                fields[ix] = tuple3_to_cell(sigma, sigma_stat, sigma_sys);
            }

            // Theory JSON for this (label, row, bin-mean kinematics).
            // We only compute theory for "10.6 GeV" and "10.2 GeV".
            auto th_it = theory_json.find(L);
            if (th_it != theory_json.end()) {
                if (xB_mid > 0.0 && Q2_mid > 0.0 && t_mid > 0.0) {
                    double Ebeam = th_it->second.value("Ebeam", (L == "10.2 GeV" ? 10.2 : 10.6));
                    json rowj;
                    rowj["row_index"] = row;
                    rowj["xB_mid"]    = xB_mid;
                    rowj["Q2_mid"]    = Q2_mid;
                    rowj["t_abs_mid"] = t_mid;

                    std::vector<double> bh_unpol(Nphi), bh_pos(Nphi), bh_neg(Nphi);
                    std::vector<double> km_unpol(Nphi), km_pos(Nphi), km_neg(Nphi);
                    std::vector<double> vgg_unpol(Nphi), vgg_pos(Nphi), vgg_neg(Nphi);

                    for (int i = 0; i < Nphi; ++i) {
                        double frac    = (i + 0.5) / static_cast<double>(Nphi);
                        double phi_rad = 2.0 * M_PI * frac;

                        // BH (helicity independent)
                        double bh_val = eval_bh_xs(Ebeam, xB_mid, Q2_mid, t_mid, phi_rad);
                        bh_unpol[i] = bh_val;
                        bh_pos[i]   = bh_val;
                        bh_neg[i]   = bh_val;

                        // KM
                        km_unpol[i] = eval_km_xs(Ebeam, xB_mid, Q2_mid, t_mid, phi_rad, "unpol");
                        km_pos[i]   = eval_km_xs(Ebeam, xB_mid, Q2_mid, t_mid, phi_rad, "pos");
                        km_neg[i]   = eval_km_xs(Ebeam, xB_mid, Q2_mid, t_mid, phi_rad, "neg");

                        // VGG
                        vgg_unpol[i] = eval_vgg_xs(Ebeam, xB_mid, Q2_mid, t_mid, phi_rad, "unpol");
                        vgg_pos[i]   = eval_vgg_xs(Ebeam, xB_mid, Q2_mid, t_mid, phi_rad, "pos");
                        vgg_neg[i]   = eval_vgg_xs(Ebeam, xB_mid, Q2_mid, t_mid, phi_rad, "neg");
                    }

                    rowj["BH"]["unpol"] = bh_unpol;
                    rowj["BH"]["pos"]   = bh_pos;
                    rowj["BH"]["neg"]   = bh_neg;

                    rowj["KM"]["unpol"] = km_unpol;
                    rowj["KM"]["pos"]   = km_pos;
                    rowj["KM"]["neg"]   = km_neg;

                    rowj["VGG"]["unpol"] = vgg_unpol;
                    rowj["VGG"]["pos"]   = vgg_pos;
                    rowj["VGG"]["neg"]   = vgg_neg;

                    std::string row_key = std::to_string(row);
                    th_it->second["rows"][row_key] = rowj;
                }
            }
        }

        out_lines.push_back(join_csv_line(fields));
    }

    // Write theory JSONs once per label
    for (const auto &L : theory_labels) {
        auto it = theory_json.find(L);
        if (it == theory_json.end()) continue;

        const std::string period_dir = canonical_period_dir(L);
        fs::path outdir = fs::path(theory_root) / period_dir;
        ensure_dir(outdir);

        fs::path fpath = outdir / "xs_phi_all.json";
        std::ofstream ofs_json(fpath);
        if (!ofs_json) {
            std::cerr << "[cross_sections] ERROR: cannot open "
                      << fpath.string() << " for writing theory JSON.\n";
        } else {
            ofs_json << std::setw(2) << it->second << "\n";
            std::cout << "[cross_sections] Wrote theory JSON for label \""
                      << L << "\" to " << fpath.string() << "\n";
        }
    }

    // Atomic write-back of CSV
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

using QTKey = std::pair<Range,Range>; // first = Q2 range, second = |t| range

struct XSGroupByXB {
    // For a given xB range, we have a map of (Q2_range, t_range) -> BinData
    std::map<QTKey, BinData> bins;
    // Optional xB index if present in CSV
    int xb_index = -1;
};

// Theory curves loaded from JSON
struct TheoryCurves {
    std::vector<double> phi_deg;
    std::vector<double> bh_unpol, bh_pos, bh_neg;
    std::vector<double> km_unpol, km_pos, km_neg;
    std::vector<double> vgg_unpol, vgg_pos, vgg_neg;
};

// Load single xs_phi_all.json for a label
static std::map<size_t, TheoryCurves>
load_theory_for_label(const std::string &label,
                      const std::string &theory_root) {
    std::map<size_t, TheoryCurves> out;

    fs::path dir  = fs::path(theory_root) / canonical_period_dir(label);
    fs::path file = dir / "xs_phi_all.json";
    if (!fs::exists(file)) {
        std::cerr << "[cross_sections] WARNING: no theory JSON file for label \""
                  << label << "\" at " << file.string() << "\n";
        return out;
    }

    std::ifstream ifs(file);
    if (!ifs) {
        std::cerr << "[cross_sections] WARNING: failed to open theory JSON for label \""
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
                  << label << "\" has no 'rows' object.\n";
        return out;
    }

    for (auto it = j["rows"].begin(); it != j["rows"].end(); ++it) {
        const std::string row_key = it.key();
        const json &cell = it.value();

        size_t row_index = 0;
        try {
            row_index = static_cast<size_t>(std::stoul(row_key));
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
              << "\" from " << file.string()
              << " with " << out.size() << " rows.\n";

    return out;
}

// Compute a sensible (ymin,ymax) for one xB slice (data + theory) for log scale
static std::pair<double,double> compute_yrange_for_slice(
    const std::map<QTKey, BinData> &bins_for_xB,
    const std::map<size_t, TheoryCurves> &theory)
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

    // Data
    for (const auto &kv : bins_for_xB) {
        const BinData &b = kv.second;
        update_from_points(b.unpol);
        update_from_points(b.pos);
        update_from_points(b.neg);

        // Theory (BH/KM/VGG, all helicities)
        if (b.have_theory_row) {
            auto it_th = theory.find(b.theory_row);
            if (it_th != theory.end()) {
                const TheoryCurves &tc = it_th->second;
                update_from_curve(tc.bh_unpol);
                update_from_curve(tc.km_unpol);
                update_from_curve(tc.km_pos);
                update_from_curve(tc.km_neg);
                update_from_curve(tc.vgg_unpol);
                update_from_curve(tc.vgg_pos);
                update_from_curve(tc.vgg_neg);
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

    // Kinematics columns
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

    // Phi: prefer phiavg,<label> if present, else phimin/phimax
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

    // Cross section columns for this label
    int c_xs_unpol = -1, c_xs_pos = -1, c_xs_neg = -1;
    try {
        c_xs_unpol = find_col(header,
            "cross sections, ep->epg, exp, " + label + ", unpol");
        c_xs_pos   = find_col(header,
            "cross sections, ep->epg, exp, " + label + ", pos");
        c_xs_neg   = find_col(header,
            "cross sections, ep->epg, exp, " + label + ", neg");
    } catch (const std::exception &e) {
        std::cerr << "[cross_sections] FATAL: " << e.what()
                  << " in plot_cross_sections_for_label.\n";
        return false;
    }

    // Build structure: xB range -> (Q2,t) -> BinData
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

        // Phi center
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
            // No cross sections for this label in this row
            continue;
        }

        XSGroupByXB &group = by_xb[xb_range];

        // Record xB index if present (first occurrence wins)
        if (group.xb_index < 0 && c_xb_idx >= 0) {
            group.xb_index = std::atoi(trim(unquote(fields[c_xb_idx])).c_str());
        }

        QTKey key(q2_range, t_range);
        BinData &bin = group.bins[key];

        // For theory lookup later
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

    // Load theory curves for this label (keyed by row index)
    std::map<size_t, TheoryCurves> theory = load_theory_for_label(label, theory_json_root);

    // ROOT styling
    gROOT->SetBatch(true);
    gStyle->SetOptStat(0);
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetLabelFont(42, "XYZ");
    gStyle->SetTextFont(42);
    gStyle->SetLineWidth(2);
    gStyle->SetFrameLineWidth(2);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    // Output dir
    fs::path outdir = fs::path(out_root_dir) / canonical_period_dir(label);
    ensure_dir(outdir);

    // Iterate xB slices in ascending xBmin order
    int xb_canvas_counter = 0;
    for (const auto &kv_xb : by_xb) {
        const Range &xb_range = kv_xb.first;
        const XSGroupByXB &group = kv_xb.second;
        const auto &bins_for_xB = group.bins;

        if (bins_for_xB.empty()) continue;

        // Collect unique Q2 and t ranges present in this xB slice
        std::set<Range> q2_set, t_set;
        for (const auto &kv : bins_for_xB) {
            const QTKey &qt = kv.first;
            q2_set.insert(qt.first);
            t_set.insert(qt.second);
        }
        std::vector<Range> q2_slice(q2_set.begin(), q2_set.end());
        std::vector<Range> t_slice(t_set.begin(), t_set.end());
        if (q2_slice.empty() || t_slice.empty()) continue;

        const int ncols = static_cast<int>(q2_slice.size());
        const int nrows = static_cast<int>(t_slice.size());

        // Per-slice y-range (data + theory), log scale
        double ymin_canvas = 1e-4;
        double ymax_canvas = 1.0;
        {
            auto yr = compute_yrange_for_slice(bins_for_xB, theory);
            ymin_canvas = yr.first;
            ymax_canvas = yr.second;
        }

        const int W = 280 * ncols + 160;
        const int H = 260 * nrows + 220;

        std::ostringstream cname;
        int xb_idx_for_name = (group.xb_index >= 0 ? group.xb_index : xb_canvas_counter);
        cname << "c_xsec_" << canonical_period_dir(label) << "_xB" << xb_idx_for_name;

        TCanvas *c = new TCanvas(cname.str().c_str(), cname.str().c_str(), W, H);

        // Top title and legend pad (slightly taller to hold legend)
        TPad *pTop = new TPad("pTop", "pTop", 0.0, 0.85, 1.0, 1.0);
        pTop->SetFillStyle(0);
        pTop->SetBorderSize(0);
        pTop->Draw();

        // Grid pad
        TPad *pGrid = new TPad("pGrid", "pGrid", 0.0, 0.00, 1.0, 0.85);
        pGrid->SetFillStyle(0);
        pGrid->SetBorderSize(0);
        pGrid->Draw();
        pGrid->cd();
        pGrid->Divide(ncols, nrows, 0.0001, 0.0001);

        // Title text + single legend on top pad
        pTop->cd();
        TLatex head;
        head.SetNDC();
        head.SetTextAlign(22);
        head.SetTextFont(42);
        head.SetTextSize(0.22);

        std::ostringstream tit;
        tit << "Cross sections, ep #rightarrow ep#gamma   " << label
            << "   x_{B} in ("
            << std::fixed << std::setprecision(3)
            << xb_range.first << ", " << xb_range.second << ")";

        head.DrawLatex(0.5, 0.70, tit.str().c_str());

        // Canvas-level legend (no per-subpad legends)
        TLegend *legTop = new TLegend(0.02, 0.05, 0.98, 0.65);
        legTop->SetBorderSize(1);
        legTop->SetLineColor(kBlack);
        legTop->SetFillColor(kWhite);
        legTop->SetFillStyle(1001);
        legTop->SetTextFont(42);
        legTop->SetTextSize(0.12);

        // Dummy objects to register styles in the top legend
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

        legTop->AddEntry(&dummy_unpol, "data unpolarized", "lep");
        legTop->AddEntry(&dummy_pos,   "data + helicity", "lep");
        legTop->AddEntry(&dummy_neg,   "data - helicity", "lep");
        legTop->AddEntry(&dummy_bh,    "BH unpolarized",  "l");
        legTop->AddEntry(&dummy_km_unpol, "KM unpolarized", "l");
        legTop->AddEntry(&dummy_km_pos,   "KM + helicity",  "l");
        legTop->AddEntry(&dummy_km_neg,   "KM - helicity",  "l");
        legTop->AddEntry(&dummy_vgg_unpol, "VGG unpolarized", "l");
        legTop->AddEntry(&dummy_vgg_pos,   "VGG + helicity",  "l");
        legTop->AddEntry(&dummy_vgg_neg,   "VGG - helicity",  "l");

        legTop->Draw();

        // Helper to draw one graph from a vector<Point>
        auto make_graph = [](const std::vector<Point> &v,
                             int mstyle, int mcolor) -> TGraphErrors* {
            if (v.empty()) return nullptr;
            int N = static_cast<int>(v.size());
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
        };

        // For each (t,Q2) combination, draw a subpad
        for (int r = 0; r < nrows; ++r) {
            const Range &t_range = t_slice[r];
            for (int cc = 0; cc < ncols; ++cc) {
                const Range &q2_range = q2_slice[cc];

                pGrid->cd(r * ncols + cc + 1);
                gPad->SetGrid(1, 1);
                gPad->SetTopMargin(0.08);
                gPad->SetBottomMargin(0.18);
                gPad->SetLeftMargin(0.16);
                gPad->SetRightMargin(0.10);
                gPad->SetLogy(true);

                QTKey key(q2_range, t_range);
                auto it_bin = bins_for_xB.find(key);

                // Always draw a frame so the canvas layout is stable
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

                // Panel label: Q^2 and |t| ranges
                TLatex lab;
                lab.SetNDC();
                lab.SetTextSize(0.070);
                lab.SetTextAlign(11);
                lab.SetTextFont(42);
                lab.DrawLatex(
                    0.15, 0.94,
                    Form("Q^{2} in (%.2f, %.2f), |t| in (%.2f, %.2f)",
                         q2_range.first, q2_range.second,
                         t_range.first,  t_range.second)
                );

                if (it_bin == bins_for_xB.end()) {
                    // No data for this (Q2, t) in this xB slice; leave empty frame.
                    continue;
                }

                BinData bin = it_bin->second; // copy so we can sort safely

                // Sort by phi for nicer curves
                auto sort_by_phi = [](std::vector<Point> &v) {
                    std::sort(v.begin(), v.end(),
                              [](const Point &a, const Point &b) {
                                  return a.phi < b.phi;
                              });
                };
                sort_by_phi(bin.unpol);
                sort_by_phi(bin.pos);
                sort_by_phi(bin.neg);

                TGraphErrors *g_unpol = make_graph(bin.unpol, 20, kBlack);
                TGraphErrors *g_pos   = make_graph(bin.pos,   24, kRed+1);
                TGraphErrors *g_neg   = make_graph(bin.neg,   25, kBlue+1);

                if (g_unpol) g_unpol->Draw("P SAME");
                if (g_pos)   g_pos->Draw("P SAME");
                if (g_neg)   g_neg->Draw("P SAME");

                // Theory curves (BH, KM, VGG; unpol and polarized for KM/VGG)
                if (bin.have_theory_row) {
                    auto it_th = theory.find(bin.theory_row);
                    if (it_th != theory.end()) {
                        const TheoryCurves &tc = it_th->second;

                        auto draw_curve = [&](const std::vector<double> &ys,
                                              int lstyle, int lcolor) {
                            if (tc.phi_deg.size() != ys.size() || ys.empty()) return;
                            int M = static_cast<int>(ys.size());
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

                        // BH (unpolarized)
                        draw_curve(tc.bh_unpol, 2, kGreen+2);

                        // KM: unpol / + / -
                        draw_curve(tc.km_unpol, 1, kMagenta+1);
                        draw_curve(tc.km_pos,   2, kMagenta+1);
                        draw_curve(tc.km_neg,   3, kMagenta+1);

                        // VGG: unpol / + / -
                        draw_curve(tc.vgg_unpol, 1, kOrange+7);
                        draw_curve(tc.vgg_pos,   2, kOrange+7);
                        draw_curve(tc.vgg_neg,   3, kOrange+7);
                    }
                }
            }
        }

        // Save canvas
        std::ostringstream fname;
        fname << "cross_sections_" << canonical_period_dir(label)
              << "_xB_" << (group.xb_index >= 0 ? group.xb_index : xb_canvas_counter)
              << ".png";
        fs::path outpath = outdir / fname.str();
        c->SaveAs(outpath.string().c_str());

        delete c;
        ++xb_canvas_counter;
    }

    std::cout << "[cross_sections] Plotted cross sections for label " << label
              << " into " << outdir.string() << "\n";

    return true;
}