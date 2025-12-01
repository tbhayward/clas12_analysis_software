// unfolded_yields_cross_check.cpp
// ------------------------------------------------------------
// Cross-check of unfolded (acceptance-corrected) yields:
//
//   Lee/pass-1 CSV column:
//     "acceptance corrected yield, ep->epg, exp"
//
//   Hayward/pass-2 CSV column (Fa18, unpolarized):
//     "acceptance corrected yield, ep->epg, exp, Fa18, unpol"
//
// The Hayward/pass-2 column is a three-tuple "(value, stat, sys)".
// We plot the central value with the stat error as TGraphErrors.
//
// Rows are compared by row index (imported binning identical).
// Output: one ROOT canvas saved to out_dir, plus console diagnostics.
// ------------------------------------------------------------

#include "unfolded_yields_cross_check.h"

#include "load_csv.h"  // your existing CSV helper

#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TAxis.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace fs = std::filesystem;

// We assume load_csv.h provides a structure like:
//
//   struct CsvTable {
//       std::vector<std::string> header;
//       std::vector< std::vector<std::string> > rows;
//   };
//
//   bool load_csv(const std::string &path, CsvTable &out);
//
// which is what the other *_cross_check modules use.

// ---------------------- helpers -----------------------------

static std::string trim(const std::string &s) {
    std::size_t first = 0;
    while (first < s.size() && std::isspace(static_cast<unsigned char>(s[first]))) {
        ++first;
    }
    if (first == s.size()) {
        return std::string();
    }
    std::size_t last = s.size() - 1;
    while (last > first && std::isspace(static_cast<unsigned char>(s[last]))) {
        --last;
    }
    return s.substr(first, last - first + 1);
}

static int find_column_index(const std::vector<std::string> &header,
                             const std::string &name) {
    for (std::size_t i = 0; i < header.size(); ++i) {
        if (header[i] == name) {
            return static_cast<int>(i);
        }
    }

    std::cerr << "[lee] FATAL (unfold): column not found: \"" << name << "\"\n";
    std::cerr << "[lee]   Available columns:\n";
    for (std::size_t i = 0; i < header.size(); ++i) {
        std::cerr << "    [" << i << "] \"" << header[i] << "\"\n";
    }
    return -1;
}

// Parse a triple "(val, stat, sys)" or "val,stat,sys".
// On success, fills val/stat/sys and returns true.
// If parsing fails, returns false and leaves outputs unchanged.
static bool parse_triple_value_stat_sys(const std::string &cell,
                                        double &val,
                                        double &stat,
                                        double &sys) {
    std::string t = cell;
    // Strip whitespace
    t.erase(std::remove_if(t.begin(), t.end(),
                           [](unsigned char c){ return std::isspace(c); }),
            t.end());

    if (t.empty()) {
        return false;
    }

    // Strip surrounding parentheses if present
    if (!t.empty() && t.front() == '(') {
        t.erase(t.begin());
    }
    if (!t.empty() && t.back() == ')') {
        t.pop_back();
    }

    std::vector<double> parts;
    std::stringstream ss(t);
    std::string token;
    while (std::getline(ss, token, ',')) {
        if (token.empty()) {
            continue;
        }
        try {
            parts.push_back(std::stod(token));
        } catch (const std::exception &) {
            return false;
        }
    }

    if (parts.size() < 2) {
        // Need at least value + stat
        return false;
    }

    val  = parts[0];
    stat = parts[1];
    sys  = (parts.size() > 2) ? parts[2] : 0.0;
    return true;
}

// Parse a "single" value (Lee column):
// - If it's already numeric, just stod.
// - If it's in triple form, take the first component.
static bool parse_single_or_first_component(const std::string &cell,
                                            double &val_out) {
    std::string t = trim(cell);
    if (t.empty()) {
        val_out = 0.0;
        return true;
    }

    // Try direct stod first
    try {
        val_out = std::stod(t);
        return true;
    } catch (const std::exception &) {
        // Fall back to triple parsing and use first component
        double v = 0.0;
        double s = 0.0;
        double y = 0.0;
        if (parse_triple_value_stat_sys(t, v, s, y)) {
            val_out = v;
            return true;
        }
    }

    return false;
}

// ---------------------- main API ----------------------------

void plot_unfolded_yields_cross_checks(const std::string &lee_csv,
                                       const std::string &hayward_csv,
                                       const std::string &out_dir) {
    std::cout << "[lee] unfolded-yields cross-check starting...\n";
    std::cout << "[lee]   Lee CSV     : " << lee_csv     << "\n";
    std::cout << "[lee]   Hayward CSV : " << hayward_csv << "\n";

    // Make sure the output directory exists
    try {
        fs::create_directories(out_dir);
    } catch (const std::exception &e) {
        std::cerr << "[lee] WARNING (unfold): failed to create output dir "
                  << out_dir << " : " << e.what() << "\n";
    }

    CsvTable lee;
    CsvTable hay;

    if (!load_csv(lee_csv, lee)) {
        std::cerr << "[lee] ERROR (unfold): failed to load Lee CSV.\n";
        return;
    }
    if (!load_csv(hayward_csv, hay)) {
        std::cerr << "[lee] ERROR (unfold): failed to load Hayward CSV.\n";
        return;
    }

    const std::size_t n_rows_lee = lee.rows.size();
    const std::size_t n_rows_hay = hay.rows.size();

    if (n_rows_lee == 0 || n_rows_hay == 0) {
        std::cerr << "[lee] FATAL (unfold): one CSV has zero rows.\n";
        std::cerr << "[lee]   Lee rows     : " << n_rows_lee << "\n";
        std::cerr << "[lee]   Hayward rows : " << n_rows_hay << "\n";
        return;
    }

    const std::size_t n_rows = std::min(n_rows_lee, n_rows_hay);
    if (n_rows_lee != n_rows_hay) {
        std::cerr << "[lee] WARNING (unfold): row counts differ (Lee="
                  << n_rows_lee << ", Hayward=" << n_rows_hay
                  << "). Using " << n_rows << " rows.\n";
    }

    const std::string lee_col_name   = "acceptance corrected yield, ep->epg, exp";
    const std::string hay_col_name   = "acceptance corrected yield, ep->epg, exp, Fa18, unpol";

    const int idx_lee = find_column_index(lee.header, lee_col_name);
    if (idx_lee < 0) {
        std::cerr << "[lee] ERROR (unfold): required Lee column missing.\n";
        return;
    }

    const int idx_hay = find_column_index(hay.header, hay_col_name);
    if (idx_hay < 0) {
        std::cerr << "[lee] ERROR (unfold): required Hayward column missing.\n";
        return;
    }

    std::vector<double> x;
    std::vector<double> lee_val;
    std::vector<double> lee_err;  // set to zero
    std::vector<double> hay_val;
    std::vector<double> hay_stat; // from tuple
    std::vector<double> hay_sys;  // unused here, but parsed for completeness

    x.reserve(n_rows);
    lee_val.reserve(n_rows);
    lee_err.reserve(n_rows);
    hay_val.reserve(n_rows);
    hay_stat.reserve(n_rows);
    hay_sys.reserve(n_rows);

    double max_y = 0.0;
    double max_diff = 0.0;
    std::size_t max_diff_row = 0;

    for (std::size_t i = 0; i < n_rows; ++i) {
        const auto &row_lee = lee.rows[i];
        const auto &row_hay = hay.rows[i];

        if (idx_lee >= static_cast<int>(row_lee.size()) ||
            idx_hay >= static_cast<int>(row_hay.size())) {
            std::cerr << "[lee] WARNING (unfold): row " << (i + 1)
                      << " has insufficient columns; skipping.\n";
            continue;
        }

        const std::string &cell_lee = row_lee[idx_lee];
        const std::string &cell_hay = row_hay[idx_hay];

        double v_lee = 0.0;
        if (!parse_single_or_first_component(cell_lee, v_lee)) {
            std::cerr << "[lee] WARNING (unfold): non-numeric Lee cell at row "
                      << (i + 1) << " : \"" << cell_lee << "\"; treating as 0.\n";
            v_lee = 0.0;
        }

        double v_hay = 0.0;
        double e_stat = 0.0;
        double e_sys  = 0.0;
        if (!parse_triple_value_stat_sys(cell_hay, v_hay, e_stat, e_sys)) {
            // Fallback: try to parse as single numeric
            if (!parse_single_or_first_component(cell_hay, v_hay)) {
                std::cerr << "[lee] WARNING (unfold): non-numeric Hayward cell at row "
                          << (i + 1) << " : \"" << cell_hay << "\"; treating as 0.\n";
                v_hay   = 0.0;
                e_stat  = 0.0;
                e_sys   = 0.0;
            } else {
                e_stat = 0.0;
                e_sys  = 0.0;
            }
        }

        const double diff = v_hay - v_lee;
        if (std::fabs(diff) > max_diff) {
            max_diff = std::fabs(diff);
            max_diff_row = i + 1; // 1-based
        }

        max_y = std::max(max_y, std::max(v_lee, v_hay));

        x.push_back(static_cast<double>(i + 1)); // bin index on x-axis
        lee_val.push_back(v_lee);
        lee_err.push_back(0.0);   // we don't have Lee stat errors here
        hay_val.push_back(v_hay);
        hay_stat.push_back(e_stat); // this is what you asked to plot
        hay_sys.push_back(e_sys);
    }

    const int n_points = static_cast<int>(x.size());
    if (n_points <= 0) {
        std::cerr << "[lee] FATAL (unfold): no valid rows to plot.\n";
        return;
    }

    std::cout << "[lee] unfolded-yields rows used: " << n_points << "\n";
    std::cout << "[lee] max |Hayward - Lee| = " << max_diff
              << " at row " << max_diff_row << "\n";

    // ----------------- Build TGraphErrors --------------------

    TGraphErrors *g_lee = new TGraphErrors(n_points);
    TGraphErrors *g_hay = new TGraphErrors(n_points);

    for (int i = 0; i < n_points; ++i) {
        g_lee->SetPoint(i, x[i], lee_val[i]);
        g_lee->SetPointError(i, 0.0, lee_err[i]);

        g_hay->SetPoint(i, x[i], hay_val[i]);
        g_hay->SetPointError(i, 0.0, hay_stat[i]);
    }

    g_lee->SetMarkerStyle(20);
    g_lee->SetMarkerSize(0.9);
    g_lee->SetMarkerColor(kBlue + 1);
    g_lee->SetLineColor(kBlue + 1);

    g_hay->SetMarkerStyle(24);
    g_hay->SetMarkerSize(0.9);
    g_hay->SetMarkerColor(kRed + 1);
    g_hay->SetLineColor(kRed + 1);

    // ----------------- Canvas and aesthetics -----------------

    gStyle->SetOptStat(0);

    TCanvas *c = new TCanvas("c_unfolded_yields",
                             "Unfolded acceptance-corrected yields cross-check",
                             1200, 800);

    c->SetGrid();

    g_hay->Draw("AP");
    g_hay->GetXaxis()->SetTitle("Bin index");
    g_hay->GetYaxis()->SetTitle("Acceptance-corrected yield");

    g_hay->GetXaxis()->CenterTitle(true);
    g_hay->GetYaxis()->CenterTitle(true);

    g_hay->GetXaxis()->SetTitleSize(0.05);
    g_hay->GetYaxis()->SetTitleSize(0.05);
    g_hay->GetXaxis()->SetLabelSize(0.045);
    g_hay->GetYaxis()->SetLabelSize(0.045);

    double ymin = 0.0;
    double ymax = (max_y <= 0.0) ? 1.0 : 1.2 * max_y;
    g_hay->GetYaxis()->SetRangeUser(ymin, ymax);

    g_lee->Draw("P SAME");

    TLegend *leg = new TLegend(0.15, 0.75, 0.55, 0.90);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.04);
    leg->AddEntry(g_lee, "Lee pass-1: acc. corr. yield (exp)", "p");
    leg->AddEntry(g_hay, "Hayward pass-2: unfolded yield (Fa18, unpol)", "pe");
    leg->Draw();

    TLatex label;
    label.SetNDC(true);
    label.SetTextSize(0.045);
    label.DrawLatex(0.15, 0.93,
                    "Unfolded acceptance-corrected yields: Lee vs Hayward (Fa18, unpol)");

    // ----------------- Save plot -----------------------------

    std::string out_png = (fs::path(out_dir) / "unfolded_yields_fa18_unpol.png").string();
    c->SaveAs(out_png.c_str());

    std::cout << "[lee] unfolded-yields cross-check plot saved to " << out_png << "\n";

    // ROOT will take ownership on exit; explicit delete is optional in this short-lived program.
}