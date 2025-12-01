#include "unfolded_yields_cross_check.h"

#include "load_csv.h"  // CsvTable + load_csv

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
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

// -------------------------
// Local helpers
// -------------------------

// Find a column index by exact header match; throw if not found
static int find_column_index(const CsvTable &table, const std::string &name) {
    for (std::size_t i = 0; i < table.header.size(); ++i) {
        if (table.header[i] == name) {
            return static_cast<int>(i);
        }
    }
    std::cerr << "[unfolded_yields] ERROR: column \"" << name
              << "\" not found in CSV header.\n";
    std::cerr << "[unfolded_yields] Header columns:\n";
    for (std::size_t i = 0; i < table.header.size(); ++i) {
        std::cerr << "  [" << i << "] " << table.header[i] << "\n";
    }
    throw std::runtime_error("Required column not found: " + name);
}

// Trim whitespace from both ends (simple)
static std::string trim(const std::string &s) {
    std::size_t i0 = 0;
    while (i0 < s.size() && std::isspace(static_cast<unsigned char>(s[i0]))) {
        ++i0;
    }
    std::size_t i1 = s.size();
    while (i1 > i0 && std::isspace(static_cast<unsigned char>(s[i1 - 1]))) {
        --i1;
    }
    return s.substr(i0, i1 - i0);
}

// Simple comma split (no quotes handling needed for our tuple)
static std::vector<std::string> split_commas(const std::string &s) {
    std::vector<std::string> out;
    std::string current;
    for (char c : s) {
        if (c == ',') {
            out.push_back(trim(current));
            current.clear();
        } else {
            current.push_back(c);
        }
    }
    if (!current.empty()) {
        out.push_back(trim(current));
    }
    return out;
}

struct Triple {
    double val;
    double stat;
    double sys;
};

// Parse a "(val, stat, sys)"-style triple; tolerate missing sys, extra spaces, etc.
// Returns true on success, false on failure.
static bool parse_triple(const std::string &raw, Triple &out) {
    std::string s = trim(raw);
    if (s.empty()) {
        return false;
    }

    // Strip leading/trailing parentheses if present
    if (!s.empty() && s.front() == '(' && s.back() == ')') {
        s = s.substr(1, s.size() - 2);
        s = trim(s);
    }

    if (s.empty()) {
        return false;
    }

    std::vector<std::string> parts = split_commas(s);
    if (parts.empty()) {
        return false;
    }

    try {
        out.val = std::stod(parts[0]);
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

// Parse a simple scalar double. Return true on success, false on failure.
static bool parse_double(const std::string &raw, double &val) {
    std::string s = trim(raw);
    if (s.empty()) {
        return false;
    }
    try {
        val = std::stod(s);
    } catch (const std::exception &e) {
        std::cerr << "[unfolded_yields] WARNING: failed to parse double from \""
                  << raw << "\" (" << e.what() << ")\n";
        return false;
    }
    return true;
}

// -------------------------
// Main cross-check function
// -------------------------

void plot_unfolded_yields_cross_checks(const std::string &lee_csv,
                                       const std::string &hayward_csv,
                                       const std::string &out_dir) {
    std::cout << "[unfolded_yields] Loading CSVs...\n";

    CsvTable lee;
    CsvTable hay;

    if (!load_csv(lee_csv, lee)) {
        std::cerr << "[unfolded_yields] ERROR: failed to load Lee CSV from "
                  << lee_csv << "\n";
        return;
    }
    if (!load_csv(hayward_csv, hay)) {
        std::cerr << "[unfolded_yields] ERROR: failed to load Hayward CSV from "
                  << hayward_csv << "\n";
        return;
    }

    // Required columns
    const std::string col_lee =
        "acceptance corrected yield, ep->epg, exp";
    const std::string col_hay =
        "acceptance corrected yield, ep->epg, exp, Fa18, unpol";

    int idx_lee = -1;
    int idx_hay = -1;

    try {
        idx_lee = find_column_index(lee, col_lee);
        idx_hay = find_column_index(hay, col_hay);
    } catch (const std::exception &ex) {
        std::cerr << "[unfolded_yields] FATAL: " << ex.what() << "\n";
        return;
    }

    const std::size_t n_rows_lee = lee.rows.size();
    const std::size_t n_rows_hay = hay.rows.size();
    const std::size_t n_rows     = std::min(n_rows_lee, n_rows_hay);

    std::cout << "[unfolded_yields] rows(Lee) = " << n_rows_lee
              << ", rows(Hayward) = " << n_rows_hay
              << ", using n_rows = " << n_rows << "\n";

    std::vector<double> x;
    std::vector<double> y_lee;
    std::vector<double> y_hay;
    std::vector<double> ey_lee;
    std::vector<double> ey_hay;

    x.reserve(n_rows);
    y_lee.reserve(n_rows);
    y_hay.reserve(n_rows);
    ey_lee.reserve(n_rows);
    ey_hay.reserve(n_rows);

    for (std::size_t i = 0; i < n_rows; ++i) {
        const auto &row_lee = lee.rows[i];
        const auto &row_hay = hay.rows[i];

        if (row_lee.size() <= static_cast<std::size_t>(idx_lee) ||
            row_hay.size() <= static_cast<std::size_t>(idx_hay)) {
            continue;
        }

        const std::string &s_lee = row_lee[idx_lee];
        const std::string &s_hay = row_hay[idx_hay];

        double lee_val = 0.0;
        if (!parse_double(s_lee, lee_val)) {
            continue;  // skip if Lee value is missing/unparseable
        }

        Triple hay_triple;
        if (!parse_triple(s_hay, hay_triple)) {
            continue;  // skip if Hayward triple is missing/unparseable
        }

        // Row index as x-value (1-based to match human row numbering)
        double xval = static_cast<double>(i + 1);

        x.push_back(xval);
        y_lee.push_back(lee_val);
        y_hay.push_back(hay_triple.val);

        // Lee CSV doesn't encode uncertainties for this column, so set 0
        ey_lee.push_back(0.0);

        // Pass-2: use the statistical component as requested
        ey_hay.push_back(hay_triple.stat);
    }

    if (x.empty()) {
        std::cerr << "[unfolded_yields] WARNING: no valid points to plot.\n";
        return;
    }

    // -------------------------
    // ROOT plotting
    // -------------------------
    gErrorIgnoreLevel = kWarning;
    gStyle->SetOptStat(0);

    TCanvas *c = new TCanvas("c_unfolded_yields",
                             "Unfolded acceptance-corrected yields cross-check",
                             900, 700);
    c->SetGrid();

    const int n_points = static_cast<int>(x.size());

    TGraphErrors *g_lee = new TGraphErrors(
        n_points,
        x.data(), y_lee.data(),
        nullptr,  // x-errors = 0
        ey_lee.data()
    );

    TGraphErrors *g_hay = new TGraphErrors(
        n_points,
        x.data(), y_hay.data(),
        nullptr,  // x-errors = 0
        ey_hay.data()
    );

    g_lee->SetName("g_lee_unfolded");
    g_hay->SetName("g_hay_unfolded");

    g_lee->SetMarkerStyle(20);
    g_lee->SetMarkerSize(0.9);
    g_lee->SetLineWidth(1);

    g_hay->SetMarkerStyle(24);
    g_hay->SetMarkerSize(0.9);
    g_hay->SetLineWidth(1);

    // Draw Hayward first so Lee can overlay if you prefer
    g_hay->SetTitle("Unfolded acceptance-corrected yields: Lee vs Hayward;Row index;Acceptance-corrected yield");
    g_hay->Draw("AP");
    g_lee->Draw("P SAME");

    // Legend
    TLegend *leg = new TLegend(0.15, 0.75, 0.45, 0.90);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(g_lee, "Lee pass-1 (unfolded)", "p");
    leg->AddEntry(g_hay, "Hayward pass-2 (Fa18, unpol)", "pe");
    leg->Draw();

    // Optionally annotate that Hayward includes stat uncertainties
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.035);
    latex.DrawLatex(0.15, 0.70,
        "Pass-2 error bars = stat uncertainty from (value, stat, sys)");

    // Adjust axis ranges a bit
    g_hay->GetXaxis()->SetLimits(0.0, static_cast<double>(n_points) + 1.0);

    double ymin = std::min(
        *std::min_element(y_lee.begin(), y_lee.end()),
        *std::min_element(y_hay.begin(), y_hay.end())
    );
    double ymax = std::max(
        *std::max_element(y_lee.begin(), y_lee.end()),
        *std::max_element(y_hay.begin(), y_hay.end())
    );

    // Add a little padding
    const double dy = (ymax - ymin) * 0.1;
    g_hay->GetYaxis()->SetRangeUser(ymin - dy, ymax + dy);

    // Save
    std::string out_file = out_dir + "/unfolded_yields_Fa18Inb.png";
    c->SaveAs(out_file.c_str());

    std::cout << "[unfolded_yields] Saved plot to " << out_file << "\n";

    // Clean up (ROOT will usually own them, but be explicit)
    delete c;
}