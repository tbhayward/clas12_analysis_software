// cross_sections.cpp
// -----------------------------------------------------------------------------
// Cross section calculator and plotting with BH / KM / VGG theory curves.
//
// This module:
//   1) Updates dvcs_pass2_analysis.csv:
//        - Fills integrated luminosity columns for all rows (constant per label).
//        - Computes cross sections for unpol/pos/neg per label where columns exist.
//        - Stores them as "(value, stat, sys)" with sys = 0.
//        - Propagates stat uncertainties from the acceptance corrected yield
//          and Frad (bin-volume stat is ignored here).
//   2) Computes theory curves at bin-mean kinematics for each row + label,
//      using 12 equally spaced phi points from 0 to 2*pi, and caches them
//      as JSON under output/jsons/cross_sections/<PeriodDir>.
//   3) Plots xB canvases with unpol/pos/neg cross sections plus BH/KM/VGG
//      theory curves as lines, one canvas per xB bin per label.
//
// This file does not add new columns to the CSV; it only updates existing ones.
// If required columns are missing, it fails with a clear error.
// -----------------------------------------------------------------------------

#include "cross_sections.h"

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

namespace fs = std::filesystem;
using json = nlohmann::json;

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
    return label;
}

// -----------------------------------------------------------------------------
// Beam energy helper (GeV)
// -----------------------------------------------------------------------------

static double beam_energy_for_label(const std::string &label) {
    // Approximate beam energies for each run period / group.
    // You can refine these if you want more precise values.
    if (label == "Fa18 Inb" ||
        label == "Fa18 Out" ||
        label == "Fa18 Inb Supp" ||
        label == "Fa18" ||
        label == "10.6 GeV") {
        return 10.604;  // fall 2018, nominal 10.6 GeV
    }

    if (label == "Sp18 Inb" ||
        label == "Sp18 Out" ||
        label == "Sp18") {
        return 10.599;  // spring 2018, nominal 10.6 GeV
    }

    if (label == "Sp19 Inb") {
        return 10.200;  // spring 2019, nominal 10.2 GeV
    }

    // Default fallback
    return 10.600;
}

// -----------------------------------------------------------------------------
// Luminosity map construction
//
// The Triple in LumiMap is used here as (value, stat, sys).
// For now we store the *total integrated charge* in nC as value and set
// stat = sys = 0. The CSV columns we fill are labeled "(nC)".
// -----------------------------------------------------------------------------

LumiMap build_lumi_map() {
    LumiMap m;

    // Per-period total charges (nC), summed over all runs with positive entries
    // in the integrated-luminosity text files.
    //
    // From the integrated charge files:
    //   rga_fa18_inb.txt   -> 2,463,154.109 nC
    //   rga_fa18_out.txt   -> 2,821,911.489 nC
    //   rga_sp18_out.txt   -> 1,228,024.094 nC
    //   rga_sp19_inb.txt   -> 1,338,346.168 nC
    //
    // Periods without an integrated-luminosity file are set to 0 here.

    m["Fa18 Inb"]      = {2463154.109, 0.0, 0.0};
    m["Fa18 Out"]      = {2821911.489, 0.0, 0.0};
    m["Fa18 Inb Supp"] = {0.0,        0.0, 0.0};   // no file provided
    m["Sp18 Inb"]      = {0.0,        0.0, 0.0};   // no file provided
    m["Sp18 Out"]      = {1228024.094, 0.0, 0.0};
    m["Sp19 Inb"]      = {1338346.168, 0.0, 0.0};

    // Combined groups (sums of periods)
    auto sum = [&](const std::vector<std::string> &keys) -> Triple {
        Triple r{0.0, 0.0, 0.0};
        for (const auto &k : keys) {
            auto it = m.find(k);
            if (it == m.end()) {
                std::cerr << "[cross_sections] ERROR: missing lumi for \""
                          << k << "\" while building combined groups.\n";
                continue;
            }
            r.value += it->second.value;
            r.stat   = std::hypot(r.stat, it->second.stat);
            r.sys    = std::hypot(r.sys, it->second.sys);
        }
        return r;
    };

    // Fa18 = Fa18 Inb + Fa18 Out (+ Fa18 Inb Supp if you ever fill it)
    m["Fa18"] = sum({"Fa18 Inb", "Fa18 Out", "Fa18 Inb Supp"});

    // Sp18 = Sp18 Inb + Sp18 Out
    m["Sp18"] = sum({"Sp18 Inb", "Sp18 Out"});

    // 10.6 GeV combined group = all 10.6-ish 2018 periods, not Sp19 (10.2 GeV)
    m["10.6 GeV"] = sum({"Fa18 Inb", "Fa18 Out", "Fa18 Inb Supp",
                         "Sp18 Inb", "Sp18 Out"});

    return m;
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

    if (s.front() == '(' && s.back() == ')') {
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
// Placeholder theory model interface
// -----------------------------------------------------------------------------

// Hook these to your actual model_predictions implementation.
// For now they return 0.0 so this compiles and runs.

static double eval_bh_xs(double /*Ebeam*/,
                         double /*xB*/, double /*Q2*/, double /*t_abs*/,
                         double /*phi_rad*/, const std::string &/*helicity*/) {
    return 0.0;
}

static double eval_km_xs(double /*Ebeam*/,
                         double /*xB*/, double /*Q2*/, double /*t_abs*/,
                         double /*phi_rad*/, const std::string &/*helicity*/) {
    return 0.0;
}

static double eval_vgg_xs(double /*Ebeam*/,
                          double /*xB*/, double /*Q2*/, double /*t_abs*/,
                          double /*phi_rad*/, const std::string &/*helicity*/) {
    return 0.0;
}

// -----------------------------------------------------------------------------
// JSON helpers
// -----------------------------------------------------------------------------

static void ensure_dir(const fs::path &p) {
    if (!fs::exists(p)) {
        fs::create_directories(p);
    }
}

// For each (label, row) we write a JSON file under:
//   output/jsons/cross_sections/<PeriodDir>/xs_phi_row<row>.json
//
// containing phi_deg and BH/KM/VGG arrays for unpol/pos/neg.
// Phi sampling: 12 equally spaced points from 0 to 2*pi.
static void write_theory_json_for_row(const std::string &label,
                                      const std::string &theory_root,
                                      size_t row_index_in_csv,
                                      double xB_mid,
                                      double Q2_mid,
                                      double t_abs_mid,
                                      double Ebeam) {
    const std::string period_dir = canonical_period_dir(label);
    fs::path outdir = fs::path(theory_root) / period_dir;
    ensure_dir(outdir);

    std::ostringstream fname;
    fname << "xs_phi_row" << row_index_in_csv << ".json";
    fs::path fpath = outdir / fname.str();

    const int Nphi = 12;
    std::vector<double> phi_deg(Nphi);
    std::vector<double> bh_unpol(Nphi), km_unpol(Nphi), vgg_unpol(Nphi);
    std::vector<double> bh_pos(Nphi),   km_pos(Nphi),   vgg_pos(Nphi);
    std::vector<double> bh_neg(Nphi),   km_neg(Nphi),   vgg_neg(Nphi);

    for (int i = 0; i < Nphi; ++i) {
        double frac = static_cast<double>(i) / Nphi;
        double phi_rad = 2.0 * M_PI * frac;
        double phi_d   = phi_rad * 180.0 / M_PI;
        phi_deg[i] = phi_d;

        bh_unpol[i]  = eval_bh_xs(Ebeam, xB_mid, Q2_mid, t_abs_mid, phi_rad, "unpol");
        km_unpol[i]  = eval_km_xs(Ebeam, xB_mid, Q2_mid, t_abs_mid, phi_rad, "unpol");
        vgg_unpol[i] = eval_vgg_xs(Ebeam, xB_mid, Q2_mid, t_abs_mid, phi_rad, "unpol");

        bh_pos[i]  = eval_bh_xs(Ebeam, xB_mid, Q2_mid, t_abs_mid, phi_rad, "pos");
        km_pos[i]  = eval_km_xs(Ebeam, xB_mid, Q2_mid, t_abs_mid, phi_rad, "pos");
        vgg_pos[i] = eval_vgg_xs(Ebeam, xB_mid, Q2_mid, t_abs_mid, phi_rad, "pos");

        bh_neg[i]  = eval_bh_xs(Ebeam, xB_mid, Q2_mid, t_abs_mid, phi_rad, "neg");
        km_neg[i]  = eval_km_xs(Ebeam, xB_mid, Q2_mid, t_abs_mid, phi_rad, "neg");
        vgg_neg[i] = eval_vgg_xs(Ebeam, xB_mid, Q2_mid, t_abs_mid, phi_rad, "neg");
    }

    json j;
    j["label"]      = label;
    j["row_index"]  = row_index_in_csv;
    j["xB_mid"]     = xB_mid;
    j["Q2_mid"]     = Q2_mid;
    j["t_abs_mid"]  = t_abs_mid;
    j["Ebeam"]      = Ebeam;
    j["phi_deg"]    = phi_deg;

    j["BH"]["unpol"] = bh_unpol;
    j["BH"]["pos"]   = bh_pos;
    j["BH"]["neg"]   = bh_neg;

    j["KM"]["unpol"] = km_unpol;
    j["KM"]["pos"]   = km_pos;
    j["KM"]["neg"]   = km_neg;

    j["VGG"]["unpol"] = vgg_unpol;
    j["VGG"]["pos"]   = vgg_pos;
    j["VGG"]["neg"]   = vgg_neg;

    std::ofstream ofs(fpath);
    ofs << std::setw(2) << j << "\n";
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

    // Header
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

    // Radiative and bin-volume columns
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
        "integrated luminosity, 10.6 GeV (nC)"
    };
    const std::string lumi_col_fa18_supp =
        "integrated luminosity, Fa18 Inb Supp (nC)";

    std::map<std::string, int> lumi_col_idx;
    // Required ones
    try {
        for (const auto &col_name : lumi_col_names_required) {
            int idx = find_col(header, col_name);
            lumi_col_idx[col_name] = idx;
        }
    } catch (const std::exception &e) {
        std::cerr << "[cross_sections] FATAL: " << e.what() << "\n";
        return false;
    }
    // Optional Fa18 Inb Supp
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
        "Fa18", "Sp18", "10.6 GeV"
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

    // Output lines
    std::vector<std::string> out_lines;
    out_lines.reserve(lines.size());
    out_lines.push_back(lines[0]);  // header

    // Theory JSON root
    const std::string theory_root = "output/jsons/cross_sections";
    ensure_dir(theory_root);

    // Process rows
    for (size_t row = 1; row < lines.size(); ++row) {
        if (lines[row].empty()) {
            out_lines.push_back(lines[row]);
            continue;
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

            auto it_lumi_col = lumi_col_idx.find(lumi_col_name);
            if (it_lumi_col != lumi_col_idx.end()) {
                fields[it_lumi_col->second] =
                    tuple3_to_cell(Lumi.value, Lumi.stat, Lumi.sys);
            }

            // Now compute cross sections for each helicity where we have columns
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

                double denom = vbin.value * Lumi.value;
                if (denom <= 0.0) {
                    std::cerr << "[cross_sections] WARNING: zero or negative denom "
                              << "(vbin * lumi) at row " << row
                              << " for label=" << L << ", helicity=" << h << ".\n";
                    continue;
                }

                // sigma = yield * Frad / (bin_volume * Lumi)
                double sigma = Y.value * frad.value / denom;

                // Stat error from Y and Frad (bin-volume stat omitted here)
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

            // Theory JSON for this (label, row, bin-mean kinematics)
            if (xB_mid > 0.0 && Q2_mid > 0.0 && t_mid > 0.0) {
                double Ebeam = beam_energy_for_label(L);
                write_theory_json_for_row(L, theory_root, row,
                                          xB_mid, Q2_mid, t_mid, Ebeam);
            }
        }

        (void)xb_idx_val; // suppress unused warning if not used later

        out_lines.push_back(join_csv_line(fields));
    }

    // Atomic write-back
    fs::path csv_path(csv_main);
    fs::path tmp_path = csv_path;
    tmp_path += ".tmp";

    std::ofstream ofs(tmp_path);
    if (!ofs) {
        std::cerr << "[cross_sections] ERROR: cannot open " << tmp_path
                  << " for writing.\n";
        return false;
    }
    for (const auto &l : out_lines) {
        ofs << l << "\n";
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

struct XSGroup {
    std::vector<Point> unpol;
    std::vector<Point> pos;
    std::vector<Point> neg;
    double xbmin = 0.0;
    double xbmax = 0.0;
    int xb_idx   = -1;
};

struct TheoryCurves {
    std::vector<double> phi_deg;
    std::vector<double> bh_unpol, km_unpol, vgg_unpol;
    std::vector<double> bh_pos,   km_pos,   vgg_pos;
    std::vector<double> bh_neg,   km_neg,   vgg_neg;
};

// Load theory JSONs for a label into a map keyed by row index
static std::map<size_t, TheoryCurves>
load_theory_for_label(const std::string &label,
                      const std::string &theory_root) {
    std::map<size_t, TheoryCurves> out;

    fs::path dir = fs::path(theory_root) / canonical_period_dir(label);
    if (!fs::exists(dir) || !fs::is_directory(dir)) {
        std::cerr << "[cross_sections] WARNING: no theory JSON dir for label \""
                  << label << "\" at " << dir.string() << "\n";
        return out;
    }

    for (const auto &entry : fs::directory_iterator(dir)) {
        if (!entry.is_regular_file()) continue;
        if (entry.path().extension() != ".json") continue;

        std::ifstream ifs(entry.path());
        if (!ifs) continue;
        json j;
        try {
            ifs >> j;
        } catch (...) {
            continue;
        }

        size_t row_index = j.value("row_index", static_cast<size_t>(0));

        TheoryCurves tc;
        tc.phi_deg   = j.value("phi_deg", std::vector<double>{});
        tc.bh_unpol  = j["BH"].value("unpol", std::vector<double>{});
        tc.bh_pos    = j["BH"].value("pos",   std::vector<double>{});
        tc.bh_neg    = j["BH"].value("neg",   std::vector<double>{});
        tc.km_unpol  = j["KM"].value("unpol", std::vector<double>{});
        tc.km_pos    = j["KM"].value("pos",   std::vector<double>{});
        tc.km_neg    = j["KM"].value("neg",   std::vector<double>{});
        tc.vgg_unpol = j["VGG"].value("unpol", std::vector<double>{});
        tc.vgg_pos   = j["VGG"].value("pos",   std::vector<double>{});
        tc.vgg_neg   = j["VGG"].value("neg",   std::vector<double>{});

        if (!tc.phi_deg.empty()) {
            out[row_index] = std::move(tc);
        }
    }
    return out;
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

    int c_xb_min = -1, c_xb_max = -1;
    int c_xb_idx = -1;
    try {
        c_xb_min = find_col(header, "xBmin");
        c_xb_max = find_col(header, "xBmax");
    } catch (const std::exception &e) {
        std::cerr << "[cross_sections] FATAL: " << e.what()
                  << " in plot_cross_sections_for_label.\n";
        return false;
    }
    c_xb_idx = find_col_optional(header, "xB index");

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

    // Group data rows by xB index (one canvas per xB bin)
    std::map<int, XSGroup> by_xb;

    for (size_t row = 1; row < lines.size(); ++row) {
        if (lines[row].empty()) continue;
        std::vector<std::string> fields = split_csv_line(lines[row]);
        if (fields.size() != header.size()) continue;

        int xb_index = -1;
        if (c_xb_idx >= 0) {
            xb_index = std::atoi(trim(unquote(fields[c_xb_idx])).c_str());
        } else {
            double xbmin = std::atof(trim(unquote(fields[c_xb_min])).c_str());
            xb_index = static_cast<int>(std::round(1000.0 * xbmin));
        }

        XSGroup &g = by_xb[xb_index];
        g.xb_idx = xb_index;
        g.xbmin = std::atof(trim(unquote(fields[c_xb_min])).c_str());
        g.xbmax = std::atof(trim(unquote(fields[c_xb_max])).c_str());

        double phi = 0.0;
        if (c_phiavg >= 0) {
            phi = std::atof(trim(unquote(fields[c_phiavg])).c_str());
        } else {
            double pmin = std::atof(trim(unquote(fields[c_phimin])).c_str());
            double pmax = std::atof(trim(unquote(fields[c_phimax])).c_str());
            phi = 0.5 * (pmin + pmax);
        }

        auto add_point = [&](int idx, std::vector<Point> &vec) {
            if (idx < 0) return;
            Triple t = parse_tuple3(fields[idx]);
            if (t.value <= 0.0) return;
            Point p;
            p.phi    = phi;
            p.xs     = t.value;
            p.xs_err = t.stat;
            vec.push_back(p);
        };

        add_point(c_xs_unpol, g.unpol);
        add_point(c_xs_pos,   g.pos);
        add_point(c_xs_neg,   g.neg);
    }

    // Load theory curves for this label (keyed by row index)
    std::map<size_t, TheoryCurves> theory = load_theory_for_label(label, theory_json_root);

    // ROOT styling
    gROOT->SetBatch(true);
    gStyle->SetOptStat(0);

    // Output dir
    fs::path outdir = fs::path(out_root_dir) / canonical_period_dir(label);
    ensure_dir(outdir);

    const char *helicity_titles[3] = {
        "unpolarized", "positive helicity", "negative helicity"
    };

    // Loop over xB groups and make one canvas per xB bin
    for (const auto &kv : by_xb) {
        int xb_index = kv.first;
        const XSGroup &g = kv.second;

        if (g.unpol.empty() && g.pos.empty() && g.neg.empty()) continue;

        // Sort by phi
        auto sort_by_phi = [](std::vector<Point> &v) {
            std::sort(v.begin(), v.end(),
                      [](const Point &a, const Point &b) { return a.phi < b.phi; });
        };
        XSGroup gg = g;  // copy to sort
        sort_by_phi(gg.unpol);
        sort_by_phi(gg.pos);
        sort_by_phi(gg.neg);

        // Canvas 1x3: unpol, pos, neg
        std::ostringstream cname;
        cname << "c_xs_" << canonical_period_dir(label) << "_xb" << xb_index;
        TCanvas *c = new TCanvas(cname.str().c_str(), "", 1600, 600);
        c->Divide(3, 1);

        const std::vector<Point> *data_vecs[3] = {
            &gg.unpol, &gg.pos, &gg.neg
        };

        // We will pick one representative TheoryCurves object to plot BH/KM/VGG.
        // Here we just use the first in the map if any exist.
        const TheoryCurves *tc_ptr = nullptr;
        if (!theory.empty()) {
            tc_ptr = &theory.begin()->second;
        }

        for (int ipad = 0; ipad < 3; ++ipad) {
            c->cd(ipad + 1);
            gPad->SetLeftMargin(0.16);
            gPad->SetRightMargin(0.06);
            gPad->SetBottomMargin(0.18);
            gPad->SetTopMargin(0.12);
            gPad->SetGridx(true);
            gPad->SetGridy(true);

            const std::vector<Point> &v = *data_vecs[ipad];
            if (v.empty()) continue;

            int N = static_cast<int>(v.size());
            std::vector<double> x(N), y(N), ex(N), ey(N);
            for (int i = 0; i < N; ++i) {
                x[i]  = v[i].phi;
                y[i]  = v[i].xs;
                ex[i] = 0.0;
                ey[i] = v[i].xs_err;
            }

            TGraphErrors *gdata = new TGraphErrors(N, x.data(), y.data(),
                                                   ex.data(), ey.data());
            gdata->SetMarkerStyle(20);
            gdata->SetMarkerSize(1.0);
            gdata->SetLineWidth(1);
            gdata->SetTitle("");
            gdata->GetXaxis()->SetTitle("#phi (deg)");
            gdata->GetYaxis()->SetTitle("d^{4}#sigma / (dx_{B} dQ^{2} d|t| d#phi)");
            gdata->GetXaxis()->CenterTitle();
            gdata->GetYaxis()->CenterTitle();
            gdata->GetXaxis()->SetTitleSize(0.07);
            gdata->GetYaxis()->SetTitleSize(0.07);
            gdata->GetXaxis()->SetLabelSize(0.06);
            gdata->GetYaxis()->SetLabelSize(0.06);
            gdata->GetXaxis()->SetNdivisions(505);
            gdata->Draw("AP");

            double ymax = 0.0;
            for (int i = 0; i < N; ++i) {
                double val = y[i] + ey[i];
                if (val > ymax) ymax = val;
            }
            if (ymax <= 0.0) ymax = 1.0;
            gdata->GetYaxis()->SetRangeUser(0.0, 1.15 * ymax);

            // Title: period label + xB range + helicity
            std::ostringstream tss;
            tss << label << ", x_{B} in ("
                << std::fixed << std::setprecision(3)
                << g.xbmin << ", " << g.xbmax << ")"
                << "  (" << helicity_titles[ipad] << ")";

            TLatex *lat = new TLatex();
            lat->SetNDC(true);
            lat->SetTextSize(0.06);
            lat->DrawLatex(0.15, 0.90, tss.str().c_str());

            TLegend *leg = new TLegend(0.55, 0.65, 0.93, 0.93);
            leg->SetFillStyle(1001);
            leg->SetFillColor(kWhite);
            leg->SetBorderSize(1);
            leg->SetTextSize(0.055);
            leg->AddEntry(gdata, "data", "pe");

            // Theory curves for this label (single representative row)
            if (tc_ptr != nullptr) {
                const TheoryCurves &tc = *tc_ptr;

                auto draw_curve = [&](const std::vector<double> &ys,
                                      int lstyle, int lcolor,
                                      const char *label_curve) {
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
                    leg->AddEntry(gth, label_curve, "l");
                };

                if (ipad == 0) {
                    draw_curve(tc.bh_unpol, 2, kRed + 1,    "BH (unpol)");
                    draw_curve(tc.km_unpol, 3, kBlue + 1,   "KM (unpol)");
                    draw_curve(tc.vgg_unpol, 4, kGreen + 2, "VGG (unpol)");
                } else if (ipad == 1) {
                    draw_curve(tc.bh_pos,   2, kRed + 1,    "BH (pos)");
                    draw_curve(tc.km_pos,   3, kBlue + 1,   "KM (pos)");
                    draw_curve(tc.vgg_pos,  4, kGreen + 2,  "VGG (pos)");
                } else {
                    draw_curve(tc.bh_neg,   2, kRed + 1,    "BH (neg)");
                    draw_curve(tc.km_neg,   3, kBlue + 1,   "KM (neg)");
                    draw_curve(tc.vgg_neg,  4, kGreen + 2,  "VGG (neg)");
                }
            }

            leg->Draw("SAME");
        }

        // Save canvas
        std::ostringstream fname;
        fname << "cross_sections_" << canonical_period_dir(label)
              << "_xB_" << xb_index << ".png";
        fs::path outpath = outdir / fname.str();
        c->SaveAs(outpath.string().c_str());

        delete c;
    }

    return true;
}