// norm_cross_sections.cpp
// -----------------------------------------------------------------------------
// Apply overall-normalization factors to already-computed DVCS cross sections
// and make the same plots as cross_sections.cpp, but using the normalized values.
//
// Convention used here:
//
//   sigma_normed = N * sigma_raw
//
// where N is read from:
//   "overall normalization, ep->epg, exp, <Label>, <hel>"
//
// and sigma_raw is read from:
//   "cross sections, ep->epg, exp, <Label>, <hel>"
//
// The normalized result is written to:
//   "normed cross sections, ep->epg, exp, <Label>, <hel>"
//
// Notes:
// - "overall normalization, ..." may be either a single number or a tuple3.
//   If it is a tuple3, we interpret it as (value, stat, sys) for N.
// - Errors are propagated assuming independence:
//
//   sigma_out = N * sigma
//   stat_out  = sqrt( (N*stat_sigma)^2 + (sigma*stat_N)^2 )
//   sys_out   = sqrt( (N*sys_sigma)^2  + (sigma*sys_N)^2  )
//
// - If either the sigma cell or N cell is blank/empty, we do not write output.
// - Fail-fast on missing required columns.
// -----------------------------------------------------------------------------

#include "norm_cross_sections.h"

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

namespace {

struct Triple {
    double value;
    double stat;
    double sys;
};

// -----------------------------------------------------------------------------
// Label helpers (match cross_sections.cpp behavior)
// -----------------------------------------------------------------------------

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

static std::string theory_energy_label_for(const std::string &label) {
    if (label == "Sp19 Inb" || label == "10.2 GeV") return "10.2 GeV";
    return "10.6 GeV";
}

static double beam_energy_for_label(const std::string &label) {
    if (label == "Sp19 Inb" || label == "10.2 GeV") return 10.2;
    return 10.6;
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

static int find_col_optional(const std::vector<std::string> &header,
                             const std::string &target) {
    for (size_t i = 0; i < header.size(); ++i) {
        if (trim(unquote(header[i])) == target) {
            return (int)i;
        }
    }
    return -1;
}

static int find_col_required(const std::vector<std::string> &header,
                             const std::string &target) {
    int idx = find_col_optional(header, target);
    if (idx < 0) {
        throw std::runtime_error("Missing required column: \"" + target + "\"");
    }
    return idx;
}

// -----------------------------------------------------------------------------
// Tuple parsing/formatting
// -----------------------------------------------------------------------------

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

    std::string s = strip_all_outer_quotes(cell);
    s = trim(s);
    if (s.empty()) return out;

    if (s.front() == '(' && s.back() == ')') {
        s = s.substr(1, s.size() - 2);
        s = trim(s);
        if (s.empty()) return out;
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
    if (!token.empty()) parts.push_back(trim(token));

    auto to_double_or_zero = [](const std::string &str) -> double {
        if (str.empty()) return 0.0;
        return std::atof(str.c_str());
    };

    if (!parts.empty())      out.value = to_double_or_zero(parts[0]);
    if (parts.size() > 1U)   out.stat  = to_double_or_zero(parts[1]);
    if (parts.size() > 2U)   out.sys   = to_double_or_zero(parts[2]);

    return out;
}

static bool cell_looks_like_tuple3(const std::string &cell) {
    std::string s = strip_all_outer_quotes(cell);
    s = trim(s);
    if (s.empty()) return false;
    if (s.front() == '(' && s.back() == ')') return true;
    if (s.find(',') != std::string::npos) return true;
    return false;
}

static bool parse_norm_cell(const std::string &cell, Triple *outN) {
    if (!outN) return false;
    *outN = Triple{0.0, 0.0, 0.0};

    std::string s = strip_all_outer_quotes(cell);
    s = trim(s);
    if (s.empty()) return false;

    // If it looks like a tuple3, parse it as (N, statN, sysN).
    if (cell_looks_like_tuple3(cell)) {
        Triple t = parse_tuple3(cell);
        if (!(t.value > 0.0)) return false;
        *outN = t;
        return true;
    }

    // Otherwise interpret as a bare scalar N.
    double N = std::atof(s.c_str());
    if (!(N > 0.0)) return false;

    outN->value = N;
    outN->stat  = 0.0;
    outN->sys   = 0.0;
    return true;
}

static std::string tuple3_to_cell(double value, double stat, double sys) {
    std::ostringstream oss;
    oss << "("
        << std::setprecision(10) << value << ", "
        << std::setprecision(10) << stat  << ", "
        << std::setprecision(10) << sys   << ")";
    return oss.str();
}

static void ensure_dir(const fs::path &p) {
    if (!fs::exists(p)) {
        fs::create_directories(p);
    }
}

// -----------------------------------------------------------------------------
// Theory loading (same structure as cross_sections.cpp)
// -----------------------------------------------------------------------------

struct TheoryCurves {
    std::vector<double> phi_deg;
    std::vector<double> bh_unpol, bh_pos, bh_neg;
    std::vector<double> km_unpol, km_pos, km_neg;
    std::vector<double> vgg_unpol, vgg_pos, vgg_neg;
};

static std::map<size_t, TheoryCurves>
load_theory_for_label(const std::string &label,
                      const std::string &theory_root) {
    std::map<size_t, TheoryCurves> out;

    std::string energy_label = theory_energy_label_for(label);
    fs::path dir  = fs::path(theory_root) / canonical_period_dir(energy_label);
    fs::path file = dir / "xs_phi_all.json";

    if (!fs::exists(file)) {
        std::cerr << "[norm_cross_sections] WARNING: no theory JSON for label \""
                  << label << "\" at " << file.string() << "\n";
        return out;
    }

    std::ifstream ifs(file);
    if (!ifs) {
        std::cerr << "[norm_cross_sections] WARNING: cannot open theory JSON for label \""
                  << label << "\" at " << file.string() << "\n";
        return out;
    }

    json j;
    try {
        ifs >> j;
    } catch (...) {
        std::cerr << "[norm_cross_sections] WARNING: malformed theory JSON for label \""
                  << label << "\" at " << file.string() << "\n";
        return out;
    }

    std::vector<double> phi_deg = j.value("phi_deg", std::vector<double>{});
    if (phi_deg.empty()) {
        std::cerr << "[norm_cross_sections] WARNING: theory JSON for label \""
                  << label << "\" has empty phi_deg.\n";
        return out;
    }

    if (!j.contains("rows") || !j["rows"].is_object()) {
        std::cerr << "[norm_cross_sections] WARNING: theory JSON for label \""
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

    std::cout << "[norm_cross_sections] Loaded theory for label \"" << label
              << "\" (energy " << energy_label << ") rows=" << out.size()
              << " from " << file.string() << "\n";

    return out;
}

// -----------------------------------------------------------------------------
// Plotting structures (copied pattern from cross_sections.cpp)
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

using QTKey = std::pair<Range, Range>;

struct XSGroupByXB {
    std::map<QTKey, BinData> bins;
    int xb_index = -1;
};

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
        if (mode == XSecPanelMode::All || mode == XSecPanelMode::UnpolOnly) update_from_points(bin->unpol);
        if (mode == XSecPanelMode::All || mode == XSecPanelMode::PosOnly)   update_from_points(bin->pos);
        if (mode == XSecPanelMode::All || mode == XSecPanelMode::NegOnly)   update_from_points(bin->neg);

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
        if (ymin <= 0.0 || !std::isfinite(ymin)) ymin = ymax * 1e-3;
        double logmin = std::pow(10.0, std::floor(std::log10(ymin)));
        double logmax = std::pow(10.0, std::ceil(std::log10(ymax)));
        ymin = std::max(1e-4, logmin);
        ymax = logmax;
    }

    return std::make_pair(ymin, 1.2 * ymax);
}

static TGraphErrors *make_xsec_graph(const std::vector<Point> &v,
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
    cname << "c_norm_xsec_";
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
    tit << "Normalized cross sections, ep #rightarrow ep#gamma   " << label
        << "   x_{B} in ("
        << std::fixed << std::setprecision(3)
        << xb_range.first << ", " << xb_range.second << ")";
    if      (mode == XSecPanelMode::UnpolOnly) tit << "   (unpolarized only)";
    else if (mode == XSecPanelMode::PosOnly)   tit << "   (+ helicity only)";
    else if (mode == XSecPanelMode::NegOnly)   tit << "   (- helicity only)";

    head.DrawLatex(0.5, 0.86, tit.str().c_str());

    TGraphErrors dummy_unpol, dummy_pos, dummy_neg;
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

    TGraph dummy_bh, dummy_km_unpol, dummy_km_pos, dummy_km_neg;
    TGraph dummy_vgg_unpol, dummy_vgg_pos, dummy_vgg_neg;

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
            if (it_bin != bins_for_xB.end()) bin_ptr = &(it_bin->second);

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
                          [](const Point &a, const Point &b) { return a.phi < b.phi; });
            };
            sort_by_phi(bin.unpol);
            sort_by_phi(bin.pos);
            sort_by_phi(bin.neg);

            if (mode == XSecPanelMode::All || mode == XSecPanelMode::UnpolOnly) {
                TGraphErrors *g_unpol = make_xsec_graph(bin.unpol, 20, kBlack);
                if (g_unpol) g_unpol->Draw("P SAME");
            }
            if (mode == XSecPanelMode::All || mode == XSecPanelMode::PosOnly) {
                TGraphErrors *g_pos = make_xsec_graph(bin.pos, 24, kRed + 1);
                if (g_pos) g_pos->Draw("P SAME");
            }
            if (mode == XSecPanelMode::All || mode == XSecPanelMode::NegOnly) {
                TGraphErrors *g_neg = make_xsec_graph(bin.neg, 25, kBlue + 1);
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
    fname << "normed_cross_sections_";
    if (mode == XSecPanelMode::UnpolOnly)      fname << "unpol_";
    else if (mode == XSecPanelMode::PosOnly)   fname << "pos_";
    else if (mode == XSecPanelMode::NegOnly)   fname << "neg_";
    fname << canonical_period_dir(label)
          << "_xB_" << xb_idx_for_name << ".png";

    fs::path outpath = outdir / fname.str();
    c->SaveAs(outpath.string().c_str());

    delete c;
}

} // end anonymous namespace

// -----------------------------------------------------------------------------
// CSV updater
// -----------------------------------------------------------------------------

bool update_normed_cross_sections_csv(const std::string &csv_main) {
    std::ifstream ifs(csv_main);
    if (!ifs) {
        std::cerr << "[norm_cross_sections] ERROR: cannot open " << csv_main
                  << " for reading.\n";
        return false;
    }

    std::vector<std::string> lines;
    std::string line;
    while (std::getline(ifs, line)) lines.push_back(line);
    ifs.close();

    if (lines.empty()) {
        std::cerr << "[norm_cross_sections] ERROR: CSV " << csv_main
                  << " is empty.\n";
        return false;
    }

    std::vector<std::string> header = split_csv_line(lines[0]);

    const std::vector<std::string> helicities = {"unpol", "pos", "neg"};
    const std::vector<std::string> labels = {
        "Fa18 Inb", "Fa18 Out",
        "Sp19 Inb",
        "Sp18 Inb", "Sp18 Out",
        "10.6 GeV",
        "Fa18",
        "Sp18"
    };

    // Required column mapping:
    // - cross sections, ...
    // - overall normalization, ...
    // - normed cross sections, ...
    struct ColTriplet {
        int xs_in  = -1;
        int norm   = -1;
        int xs_out = -1;
    };
    std::map<std::string, ColTriplet> cmap;

    std::vector<std::string> missing;

    for (const auto &L : labels) {
        for (const auto &h : helicities) {
            std::string xs_in_name =
                "cross sections, ep->epg, exp, " + L + ", " + h;
            std::string norm_name =
                "overall normalization, ep->epg, exp, " + L + ", " + h;
            std::string xs_out_name =
                "normed cross sections, ep->epg, exp, " + L + ", " + h;

            int i_xs_in  = find_col_optional(header, xs_in_name);
            int i_norm   = find_col_optional(header, norm_name);
            int i_xs_out = find_col_optional(header, xs_out_name);

            // We require all three to exist as columns (fail fast if not).
            if (i_xs_in < 0)  missing.push_back(xs_in_name);
            if (i_norm < 0)   missing.push_back(norm_name);
            if (i_xs_out < 0) missing.push_back(xs_out_name);

            if (i_xs_in >= 0 && i_norm >= 0 && i_xs_out >= 0) {
                std::string key = L + "|" + h;
                cmap[key] = ColTriplet{i_xs_in, i_norm, i_xs_out};
            }
        }
    }

    if (!missing.empty()) {
        std::ostringstream oss;
        oss << "[norm_cross_sections] FATAL: missing required columns (" << missing.size() << "):\n";
        for (const auto &m : missing) oss << "  - " << m << "\n";
        std::cerr << oss.str();
        return false;
    }

    std::vector<std::string> out_lines;
    out_lines.reserve(lines.size());
    out_lines.push_back(lines[0]);

    const size_t n_data_rows = (lines.size() > 0 ? lines.size() - 1 : 0);
    std::cout << "[norm_cross_sections] update_normed_cross_sections_csv: data rows = "
              << n_data_rows << "\n";

    int next_pct = 1;

    for (size_t row = 1; row < lines.size(); ++row) {
        if (lines[row].empty()) {
            out_lines.push_back(lines[row]);
            continue;
        }

        if (n_data_rows > 0 && next_pct <= 100) {
            double frac = 100.0 * (double)row / (double)n_data_rows;
            if (frac >= next_pct) {
                std::cout << "[norm_cross_sections] ~"
                          << next_pct << "% rows processed (row "
                          << row << " / " << n_data_rows << ")\n";
                next_pct += 10;
            }
        }

        std::vector<std::string> fields = split_csv_line(lines[row]);
        if (fields.size() != header.size()) {
            std::cerr << "[norm_cross_sections] WARNING: row " << row
                      << " has wrong number of fields, copying unchanged.\n";
            out_lines.push_back(lines[row]);
            continue;
        }

        for (const auto &L : labels) {
            for (const auto &h : helicities) {
                std::string key = L + "|" + h;
                auto it = cmap.find(key);
                if (it == cmap.end()) continue;

                const ColTriplet &C = it->second;

                // If sigma cell is blank -> do not write output.
                std::string xs_cell = trim(unquote(fields[C.xs_in]));
                if (strip_all_outer_quotes(xs_cell).empty()) {
                    continue;
                }

                Triple sigma_in = parse_tuple3(fields[C.xs_in]);
                if (!(sigma_in.value > 0.0)) {
                    continue;
                }

                // If norm cell is blank -> do not write output.
                Triple N{0.0, 0.0, 0.0};
                if (!parse_norm_cell(fields[C.norm], &N)) {
                    continue;
                }

                // Convention: sigma_out = N * sigma_in
                double sigma_out = N.value * sigma_in.value;

                double stat2 = 0.0;
                double sys2  = 0.0;

                if (sigma_in.stat > 0.0) stat2 += (N.value * sigma_in.stat) * (N.value * sigma_in.stat);
                if (N.stat > 0.0)        stat2 += (sigma_in.value * N.stat) * (sigma_in.value * N.stat);

                if (sigma_in.sys > 0.0)  sys2  += (N.value * sigma_in.sys) * (N.value * sigma_in.sys);
                if (N.sys > 0.0)         sys2  += (sigma_in.value * N.sys) * (sigma_in.value * N.sys);

                double stat_out = (stat2 > 0.0) ? std::sqrt(stat2) : 0.0;
                double sys_out  = (sys2  > 0.0) ? std::sqrt(sys2)  : 0.0;

                fields[C.xs_out] = tuple3_to_cell(sigma_out, stat_out, sys_out);
            }
        }

        out_lines.push_back(join_csv_line(fields));
    }

    fs::path csv_path(csv_main);
    fs::path tmp_path = csv_path;
    tmp_path += ".tmp";

    std::ofstream ofs(tmp_path);
    if (!ofs) {
        std::cerr << "[norm_cross_sections] ERROR: cannot open " << tmp_path
                  << " for writing.\n";
        return false;
    }
    for (const auto &lout : out_lines) ofs << lout << "\n";
    ofs.close();

    fs::rename(tmp_path, csv_path);
    std::cout << "[norm_cross_sections] Updated CSV with normed cross sections: "
              << csv_main << "\n";

    return true;
}

// -----------------------------------------------------------------------------
// Plotting entry point (normed columns)
// -----------------------------------------------------------------------------

bool plot_normed_cross_sections_for_label(const std::string &csv_main,
                                          const std::string &label,
                                          const std::string &theory_json_root,
                                          const std::string &out_root_dir) {
    std::ifstream ifs(csv_main);
    if (!ifs) {
        std::cerr << "[norm_cross_sections] ERROR: cannot open " << csv_main
                  << " for plotting.\n";
        return false;
    }

    std::vector<std::string> lines;
    std::string line;
    while (std::getline(ifs, line)) lines.push_back(line);
    ifs.close();

    if (lines.empty()) {
        std::cerr << "[norm_cross_sections] ERROR: CSV " << csv_main
                  << " is empty in plot_normed_cross_sections_for_label.\n";
        return false;
    }

    std::vector<std::string> header = split_csv_line(lines[0]);

    int c_xb_min = -1, c_xb_max = -1;
    int c_q2_min = -1, c_q2_max = -1;
    int c_t_min  = -1, c_t_max  = -1;
    try {
        c_xb_min = find_col_required(header, "xBmin");
        c_xb_max = find_col_required(header, "xBmax");
        c_q2_min = find_col_required(header, "Q2min");
        c_q2_max = find_col_required(header, "Q2max");
        c_t_min  = find_col_required(header, "t_abs_min");
        c_t_max  = find_col_required(header, "t_abs_max");
    } catch (const std::exception &e) {
        std::cerr << "[norm_cross_sections] FATAL: " << e.what()
                  << " in plot_normed_cross_sections_for_label.\n";
        return false;
    }

    int c_xb_idx = find_col_optional(header, "xB index");

    int c_phiavg = find_col_optional(header, "phiavg, " + label);
    int c_phimin = -1, c_phimax = -1;
    if (c_phiavg < 0) {
        c_phimin = find_col_optional(header, "phimin");
        c_phimax = find_col_optional(header, "phimax");
        if (c_phimin < 0 || c_phimax < 0) {
            std::cerr << "[norm_cross_sections] FATAL: no phiavg or phimin/phimax "
                      << "available for label " << label << ".\n";
            return false;
        }
    }

    int c_xs_unpol = find_col_optional(
        header, "normed cross sections, ep->epg, exp, " + label + ", unpol");
    int c_xs_pos   = find_col_optional(
        header, "normed cross sections, ep->epg, exp, " + label + ", pos");
    int c_xs_neg   = find_col_optional(
        header, "normed cross sections, ep->epg, exp, " + label + ", neg");

    if (c_xs_unpol < 0 && c_xs_pos < 0 && c_xs_neg < 0) {
        std::cerr << "[norm_cross_sections] INFO: no normed cross section columns for label "
                  << label << "; nothing to plot.\n";
        return true;
    }
    if (c_xs_unpol < 0 || c_xs_pos < 0 || c_xs_neg < 0) {
        std::cerr << "[norm_cross_sections] FATAL: incomplete normed cross section columns "
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
        Triple xs_pos   = parse_tuple3(fields[c_xs_pos]);
        Triple xs_neg   = parse_tuple3(fields[c_xs_neg]);
        if (xs_unpol.value <= 0.0 && xs_pos.value <= 0.0 && xs_neg.value <= 0.0) continue;

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
        add_point(xs_pos,   bin.pos);
        add_point(xs_neg,   bin.neg);
    }

    if (by_xb.empty()) {
        std::cerr << "[norm_cross_sections] WARNING: no normed xsec data found for label "
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

        std::set<Range> q2_set, t_set;
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

    std::cout << "[norm_cross_sections] Plotted normed cross sections for label "
              << label << " into " << outdir.string() << "\n";
    return true;
}