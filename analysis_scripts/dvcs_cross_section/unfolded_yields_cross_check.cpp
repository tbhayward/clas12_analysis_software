// unfolded_yields_cross_check.cpp
// -----------------------------------------------------------------------------
// Cross-check of unfolded / acceptance-corrected yields (Lee pass-1 vs
// Hayward pass-2) using *only* CSVs.
//
// Lee CSV (imports/all_bin_v3.csv):
//   - scalar column:
//       "acceptance corrected yield, ep->epg, exp"
//
// Hayward CSV (output/csvs/dvcs_pass2_analysis.csv):
//   - tuple column (value, stat, sys):
//       "acceptance corrected yield, ep->epg, exp, Fa18, unpol"
//
// We match rows by "bin index" (with "valid bin" == 1 in both CSVs) and build
// TGraphErrors of:
//
//   y_Lee(bin)      = Lee acceptance-corrected yield
//   y_Hayward(bin)  = Hayward acceptance-corrected yield
//   ey_Lee(bin)     = 0
//   ey_Hayward(bin) = stat component from the (value, stat, sys) tuple
//
// Output:
//   <out_dir>/unfolded_yields_Fa18.png
// -----------------------------------------------------------------------------

#include "unfolded_yields_cross_check.h"

#include <TCanvas.h>
#include <TError.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TLatex.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

// -------------------------
// Small utilities
// -------------------------

static inline void info(const std::string &s) {
    std::cout << "[unfolded_yields] " << s << std::endl;
}

static inline void fatal(const std::string &s) {
    std::cerr << "[unfolded_yields][FATAL] " << s << std::endl;
    throw std::runtime_error(s);
}

// Trim whitespace from both ends
static inline std::string trim(const std::string &s) {
    std::size_t i0 = 0;
    while (i0 < s.size() && std::isspace((unsigned char)s[i0])) {
        ++i0;
    }
    std::size_t i1 = s.size();
    while (i1 > i0 && std::isspace((unsigned char)s[i1 - 1])) {
        --i1;
    }
    return s.substr(i0, i1 - i0);
}

// -------------------------
// CSV helpers (same style as raw_signal_cross_check.cpp)
// -------------------------

// Split a CSV line into cells, handling double quotes.
static std::vector<std::string> split_csv_line(const std::string &line) {
    std::vector<std::string> out;
    std::string cur;
    cur.reserve(line.size());
    bool in_quotes = false;

    for (char c : line) {
        if (c == '"') {
            in_quotes = !in_quotes;
            continue;
        }
        if (c == ',' && !in_quotes) {
            out.push_back(cur);
            cur.clear();
        } else {
            cur.push_back(c);
        }
    }
    out.push_back(cur);
    return out;
}

static std::unordered_map<std::string,int>
build_header_index(const std::vector<std::string> &header) {
    std::unordered_map<std::string,int> m;
    for (int i = 0; i < (int)header.size(); ++i) {
        m[header[i]] = i;
    }
    return m;
}

static const std::string &get_col_ref(const std::vector<std::string> &row,
                                      const std::unordered_map<std::string,int> &idx,
                                      const std::string &name) {
    auto it = idx.find(name);
    static const std::string empty;
    if (it == idx.end()) {
        return empty;
    }
    int j = it->second;
    if (j < 0 || j >= (int)row.size()) {
        return empty;
    }
    return row[j];
}

static double ToDouble(const std::string &s) {
    std::string t = trim(s);
    if (t.empty()) return 0.0;
    char *endp = nullptr;
    double v = std::strtod(t.c_str(), &endp);
    if (endp == t.c_str()) return 0.0;
    return v;
}

static int ToInt(const std::string &s) {
    std::string t = trim(s);
    if (t.empty()) return 0;
    char *endp = nullptr;
    long v = std::strtol(t.c_str(), &endp, 10);
    if (endp == t.c_str()) return 0;
    return (int)v;
}

static void require_columns(const std::unordered_map<std::string,int> &idx,
                            const std::vector<std::string> &cols,
                            const std::string &which_csv) {
    std::vector<std::string> missing;
    for (const auto &c : cols) {
        if (idx.find(c) == idx.end()) {
            missing.push_back(c);
        }
    }
    if (!missing.empty()) {
        std::ostringstream oss;
        oss << "Missing columns in " << which_csv << ": ";
        for (std::size_t i = 0; i < missing.size(); ++i) {
            if (i > 0) oss << ", ";
            oss << '"' << missing[i] << '"';
        }
        fatal(oss.str());
    }
}

// -------------------------
// Tuple parser for Hayward triple "(val, stat, sys)"
// -------------------------

struct Triple {
    double val;
    double stat;
    double sys;
};

// Simple comma split (no quotes handling needed inside tuple)
static std::vector<std::string> split_commas(const std::string &s) {
    std::vector<std::string> out;
    std::string cur;
    for (char c : s) {
        if (c == ',') {
            out.push_back(trim(cur));
            cur.clear();
        } else {
            cur.push_back(c);
        }
    }
    if (!cur.empty()) {
        out.push_back(trim(cur));
    }
    return out;
}

static bool parse_triple(const std::string &raw, Triple &out) {
    std::string s = trim(raw);
    if (s.empty()) return false;

    // Strip surrounding parentheses if present
    if (!s.empty() && s.front() == '(' && s.back() == ')') {
        s = s.substr(1, s.size() - 2);
        s = trim(s);
    }
    if (s.empty()) return false;

    std::vector<std::string> parts = split_commas(s);
    if (parts.empty()) return false;

    try {
        out.val  = std::stod(parts[0]);
        out.stat = 0.0;
        out.sys  = 0.0;
        if (parts.size() >= 2) {
            out.stat = std::stod(parts[1]);
        }
        if (parts.size() >= 3) {
            out.sys = std::stod(parts[2]);
        }
    } catch (const std::exception &e) {
        std::cerr << "[unfolded_yields] WARNING: failed to parse triple from \""
                  << raw << "\" (" << e.what() << ")\n";
        return false;
    }
    return true;
}

// -------------------------
// Main cross-check function
// -------------------------

void plot_unfolded_yields_cross_checks(const std::string &lee_csv_path,
                                       const std::string &hayward_csv_path,
                                       const std::string &out_dir) {
    info("Loading Lee CSV: " + lee_csv_path);
    info("Loading Hayward CSV: " + hayward_csv_path);

    // --------- 1. Load Lee CSV and collect scalar yields by bin index ---------
    std::ifstream fin_lee(lee_csv_path);
    if (!fin_lee.is_open()) {
        fatal("Cannot open Lee CSV: " + lee_csv_path);
    }

    std::string header_line;
    if (!std::getline(fin_lee, header_line)) {
        fatal("Lee CSV appears empty: " + lee_csv_path);
    }
    std::vector<std::string> header_lee = split_csv_line(header_line);
    auto H_lee = build_header_index(header_lee);

    const std::string col_bin_lee  = "bin index";
    const std::string col_valid_lee = "valid bin";
    const std::string col_y_lee =
        "acceptance corrected yield, ep->epg, exp";

    require_columns(H_lee,
                    { col_bin_lee, col_valid_lee, col_y_lee },
                    "Lee CSV");

    std::map<int,double> lee_yield_by_bin;

    std::string line;
    int lee_rows_read  = 0;
    int lee_rows_kept  = 0;

    while (std::getline(fin_lee, line)) {
        ++lee_rows_read;
        if (line.empty()) continue;

        std::vector<std::string> cols = split_csv_line(line);
        const std::string &s_valid = get_col_ref(cols, H_lee, col_valid_lee);
        int valid = ToInt(s_valid);
        if (valid != 1) continue;

        int bin = ToInt(get_col_ref(cols, H_lee, col_bin_lee));
        const std::string &s_y = get_col_ref(cols, H_lee, col_y_lee);
        double y = ToDouble(s_y);

        lee_yield_by_bin[bin] = y;
        ++lee_rows_kept;
    }

    info("Lee rows read: " + std::to_string(lee_rows_read));
    info("Lee valid rows kept: " + std::to_string(lee_rows_kept));

    // --------- 2. Load Hayward CSV and collect triple yields by bin index -----
    std::ifstream fin_hay(hayward_csv_path);
    if (!fin_hay.is_open()) {
        fatal("Cannot open Hayward CSV: " + hayward_csv_path);
    }

    if (!std::getline(fin_hay, header_line)) {
        fatal("Hayward CSV appears empty: " + hayward_csv_path);
    }
    std::vector<std::string> header_hay = split_csv_line(header_line);
    auto H_hay = build_header_index(header_hay);

    const std::string col_bin_hay   = "bin index";
    const std::string col_valid_hay = "valid bin";
    const std::string col_y_hay =
        "acceptance corrected yield, ep->epg, exp, Fa18, unpol";

    require_columns(H_hay,
                    { col_bin_hay, col_valid_hay, col_y_hay },
                    "Hayward CSV");

    struct HayVal {
        double val;
        double stat;
    };

    std::map<int,HayVal> hay_yield_by_bin;

    int hay_rows_read = 0;
    int hay_rows_kept = 0;

    while (std::getline(fin_hay, line)) {
        ++hay_rows_read;
        if (line.empty()) continue;

        std::vector<std::string> cols = split_csv_line(line);
        const std::string &s_valid = get_col_ref(cols, H_hay, col_valid_hay);
        int valid = ToInt(s_valid);
        if (valid != 1) continue;

        int bin = ToInt(get_col_ref(cols, H_hay, col_bin_hay));
        const std::string &s_y = get_col_ref(cols, H_hay, col_y_hay);

        Triple t;
        if (!parse_triple(s_y, t)) {
            // If we cannot parse the triple, skip this bin.
            std::cerr << "[unfolded_yields] WARNING: skipping bin " << bin
                      << " (cannot parse tuple \"" << s_y << "\")\n";
            continue;
        }

        HayVal hv;
        hv.val  = t.val;
        hv.stat = t.stat;
        hay_yield_by_bin[bin] = hv;
        ++hay_rows_kept;
    }

    info("Hayward rows read: " + std::to_string(hay_rows_read));
    info("Hayward valid rows kept (parsed triple): " + std::to_string(hay_rows_kept));

    // --------- 3. Build intersection of bins and vectors for plotting ---------
    std::vector<int> bins;
    bins.reserve(std::min(lee_yield_by_bin.size(), hay_yield_by_bin.size()));

    for (const auto &kv : lee_yield_by_bin) {
        int bin = kv.first;
        if (hay_yield_by_bin.find(bin) != hay_yield_by_bin.end()) {
            bins.push_back(bin);
        }
    }

    if (bins.empty()) {
        std::cerr << "[unfolded_yields] WARNING: no overlapping bins between Lee and Hayward.\n";
        return;
    }

    std::sort(bins.begin(), bins.end());

    std::vector<double> x;
    std::vector<double> y_lee;
    std::vector<double> y_hay;
    std::vector<double> ey_lee;
    std::vector<double> ey_hay;

    x.reserve(bins.size());
    y_lee.reserve(bins.size());
    y_hay.reserve(bins.size());
    ey_lee.reserve(bins.size());
    ey_hay.reserve(bins.size());

    for (int bin : bins) {
        double lee_val = lee_yield_by_bin[bin];
        HayVal hv      = hay_yield_by_bin[bin];

        x.push_back((double)bin);
        y_lee.push_back(lee_val);
        y_hay.push_back(hv.val);
        ey_lee.push_back(0.0);        // Lee CSV has no stat error for this column
        ey_hay.push_back(hv.stat);    // stat uncertainty from (value, stat, sys)
    }

    info("Number of bins used for plot: " + std::to_string(x.size()));

    if (x.empty()) {
        std::cerr << "[unfolded_yields] WARNING: nothing to plot.\n";
        return;
    }

    // --------- 4. ROOT plotting ---------
    gErrorIgnoreLevel = kWarning;
    gStyle->SetOptStat(0);

    TCanvas *c = new TCanvas("c_unfolded_yields",
                             "Unfolded acceptance-corrected yields cross-check",
                             900, 700);
    c->SetGrid(1, 1);

    const int n_points = (int)x.size();

    TGraphErrors *g_lee = new TGraphErrors(
        n_points,
        x.data(), y_lee.data(),
        nullptr,          // x-errors = 0
        ey_lee.data()
    );

    TGraphErrors *g_hay = new TGraphErrors(
        n_points,
        x.data(), y_hay.data(),
        nullptr,          // x-errors = 0
        ey_hay.data()
    );

    g_lee->SetName("g_lee_unfolded");
    g_hay->SetName("g_hay_unfolded");

    // Styles: Lee = orange open, Hayward = black solid with error bars
    g_lee->SetMarkerStyle(24);
    g_lee->SetMarkerSize(0.9);
    g_lee->SetLineWidth(1);

    g_hay->SetMarkerStyle(20);
    g_hay->SetMarkerSize(0.9);
    g_hay->SetLineWidth(1);

    g_hay->SetTitle("Unfolded acceptance-corrected yields: Lee vs Hayward;Bin index;Acceptance-corrected yield");
    g_hay->Draw("AP");
    g_lee->Draw("P SAME");

    // Legend
    TLegend *leg = new TLegend(0.15, 0.75, 0.45, 0.90);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(g_hay, "Hayward pass-2 (Fa18, unpol)", "pe");
    leg->AddEntry(g_lee, "Lee pass-1 (unfolded)",        "p");
    leg->Draw();

    // Explain error bars
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.035);
    latex.DrawLatex(0.15, 0.70,
        "Pass-2 error bars = stat component from (value, stat, sys)");

    // Axis ranges
    g_hay->GetXaxis()->SetLimits(
        std::min_element(x.begin(), x.end())[0] - 1.0,
        std::max_element(x.begin(), x.end())[0] + 1.0
    );

    double ymin = std::min(
        *std::min_element(y_lee.begin(), y_lee.end()),
        *std::min_element(y_hay.begin(), y_hay.end())
    );
    double ymax = std::max(
        *std::max_element(y_lee.begin(), y_lee.end()),
        *std::max_element(y_hay.begin(), y_hay.end())
    );
    double dy = (ymax - ymin) * 0.10;
    if (dy <= 0.0) dy = std::abs(ymax) * 0.10;

    g_hay->GetYaxis()->SetRangeUser(ymin - dy, ymax + dy);

    std::string out_file = out_dir + "/unfolded_yields_Fa18.png";
    c->SaveAs(out_file.c_str());
    info("Saved plot to " + out_file);

    delete leg;
    delete c;
}