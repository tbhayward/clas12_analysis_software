#include "overall_normalization_study.h"

#include "model_predictions.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// ROOT
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TGraph2D.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TColor.h>
#include <TROOT.h>
#include <TF1.h>
#include <TFitResultPtr.h>
#include <TFitResult.h>
#include <TMatrixDSym.h>
#include <TMath.h>

#include <filesystem>

namespace {

// -----------------------------------------------------------------------------
// User-controlled constants for the new grid fits (explicit, no hidden fallback)
// -----------------------------------------------------------------------------
static const double kPhi0Deg = 0.001;

// Fit model:
//   xs(phi) = A * (1 + B cos(phi) + C cos(2 phi) + D cos(3 phi))
// with phi in degrees.
static double xs_harmonic_model_deg(double *x, double *p) {
    const double phi_deg = x[0];
    const double phi_rad = phi_deg * TMath::DegToRad();

    const double A = p[0];
    const double B = p[1];
    const double C = p[2];
    const double D = p[3];

    const double c1 = std::cos(1.0 * phi_rad);
    const double c2 = std::cos(2.0 * phi_rad);
    const double c3 = std::cos(3.0 * phi_rad);

    return A * (1.0 + B * c1 + C * c2 + D * c3);
}

struct TripleCell {
    double value;
    double stat;
    double sys;
};

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

static int find_col_optional(const std::vector<std::string> &header,
                             const std::string &target) {
    for (size_t i = 0; i < header.size(); ++i) {
        if (trim(unquote(header[i])) == target) {
            return (int)i;
        }
    }
    return -1;
}

static int require_col(const std::vector<std::string> &header,
                       const std::string &target) {
    const int idx = find_col_optional(header, target);
    if (idx < 0) {
        throw std::runtime_error("Missing required column: \"" + target + "\"");
    }
    return idx;
}

static std::string strip_all_outer_quotes(std::string s) {
    s = unquote(s);
    s = trim(s);

    bool changed = true;
    while (changed && s.size() >= 2) {
        changed = false;
        const char first = s.front();
        const char last  = s.back();
        if ((first == '"' && last == '"') || (first == '\'' && last == '\'')) {
            s = s.substr(1, s.size() - 2);
            s = trim(s);
            changed = true;
        }
    }
    return s;
}

static TripleCell parse_tuple3(const std::string &cell) {
    TripleCell out;
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
    parts.push_back(trim(token));

    auto to_double_or_zero = [](const std::string &str) -> double {
        if (str.empty()) return 0.0;
        return std::atof(str.c_str());
    };

    if (!parts.empty()) out.value = to_double_or_zero(parts[0]);
    if (parts.size() > 1U) out.stat = to_double_or_zero(parts[1]);
    if (parts.size() > 2U) out.sys  = to_double_or_zero(parts[2]);

    return out;
}

static bool finite_pos(double x) {
    return std::isfinite(x) && x > 0.0;
}

static Helicity helicity_from_string(const std::string &h) {
    if (h == "pos") return Helicity::Plus;
    if (h == "neg") return Helicity::Minus;
    return Helicity::Unpol;
}

static double beam_energy_for_label(const std::string &label) {
    if (label == "Sp19 Inb" || label == "10.2 GeV") {
        return 10.2;
    }
    return 10.6;
}

struct BestPhiRow {
    double xb_c;
    double q2_c;
    double t_c;   // positive |t|
    double phi_deg;
    double dist_to_edge;
    double xs;
    double xs_stat;
    size_t csv_row_index;

    BestPhiRow() :
        xb_c(std::numeric_limits<double>::quiet_NaN()),
        q2_c(std::numeric_limits<double>::quiet_NaN()),
        t_c(std::numeric_limits<double>::quiet_NaN()),
        phi_deg(std::numeric_limits<double>::quiet_NaN()),
        dist_to_edge(std::numeric_limits<double>::infinity()),
        xs(0.0),
        xs_stat(0.0),
        csv_row_index(0) {}
};

static double min_dist_to_0_or_360(double phi_deg) {
    const double d0   = std::fabs(phi_deg - 0.0);
    const double d360 = std::fabs(phi_deg - 360.0);
    return std::min(d0, d360);
}

static std::string cell_key_for_kin_bin(const std::string &xbmin_s,
                                        const std::string &xbmax_s,
                                        const std::string &q2min_s,
                                        const std::string &q2max_s,
                                        const std::string &tmin_s,
                                        const std::string &tmax_s) {
    std::string k;
    k.reserve(xbmin_s.size() + xbmax_s.size() + q2min_s.size() + q2max_s.size() +
              tmin_s.size() + tmax_s.size() + 16);
    k += xbmin_s; k += "|";
    k += xbmax_s; k += "|";
    k += q2min_s; k += "|";
    k += q2max_s; k += "|";
    k += tmin_s;  k += "|";
    k += tmax_s;
    return k;
}

static std::string sanitize_for_filename(const std::string &s) {
    std::string out;
    out.reserve(s.size());
    for (size_t i = 0; i < s.size(); ++i) {
        const unsigned char c = (unsigned char)s[i];
        if (std::isalnum(c)) {
            out.push_back((char)c);
        } else if (c == ' ' || c == '-' || c == '.') {
            out.push_back('_');
        } else {
            // drop other punctuation deterministically
        }
    }
    std::string collapsed;
    collapsed.reserve(out.size());
    bool prev_us = false;
    for (size_t i = 0; i < out.size(); ++i) {
        if (out[i] == '_') {
            if (!prev_us) collapsed.push_back('_');
            prev_us = true;
        } else {
            collapsed.push_back(out[i]);
            prev_us = false;
        }
    }
    if (!collapsed.empty() && collapsed.back() == '_') collapsed.pop_back();
    return collapsed;
}

static void ensure_output_dir_or_throw(const std::string &dir) {
    std::error_code ec;
    if (std::filesystem::exists(dir, ec)) {
        if (ec) {
            throw std::runtime_error("Failed to stat output directory: " + dir);
        }
        return;
    }
    if (!std::filesystem::create_directories(dir, ec)) {
        if (ec) {
            throw std::runtime_error("Failed to create output directory: " + dir);
        }
    }
}

// -----------------------------------------------------------------------------
// Legacy 1D/2D grouping by BH/KM15 (unchanged behavior)
// -----------------------------------------------------------------------------
enum class BhKmGroup {
    G_95_105 = 0,
    G_90_95_OR_105_110,
    G_85_90_OR_110_115,
    G_OUTSIDE_15PCT,
    G_INVALID
};

static BhKmGroup categorize_bh_over_km15(double r) {
    if (!std::isfinite(r) || r <= 0.0) return BhKmGroup::G_INVALID;

    if (r >= 0.95 && r <= 1.05) return BhKmGroup::G_95_105;

    if ((r >= 0.90 && r < 0.95) || (r > 1.05 && r <= 1.10)) return BhKmGroup::G_90_95_OR_105_110;

    if ((r > 0.85 && r < 0.90) || (r > 1.10 && r < 1.15)) return BhKmGroup::G_85_90_OR_110_115;

    if (r <= 0.85 || r >= 1.15) return BhKmGroup::G_OUTSIDE_15PCT;

    return BhKmGroup::G_INVALID;
}

struct GroupStyle {
    BhKmGroup group;
    int marker_style;
    int marker_color;
    std::string label;
};

static std::vector<GroupStyle> make_group_styles() {
    std::vector<GroupStyle> s;
    s.push_back({BhKmGroup::G_95_105,             20, kBlack,     "BH/KM15 in [0.95, 1.05]"});
    s.push_back({BhKmGroup::G_90_95_OR_105_110,   21, kBlue+1,    "BH/KM15 in [0.90, 0.95) or (1.05, 1.10]"});
    s.push_back({BhKmGroup::G_85_90_OR_110_115,   22, kRed+1,     "BH/KM15 in (0.85, 0.90) or (1.10, 1.15)"});
    s.push_back({BhKmGroup::G_OUTSIDE_15PCT,      24, kMagenta+1, "BH/KM15 <= 0.85 or >= 1.15"});
    return s;
}

static int style_index_for_group(BhKmGroup g, const std::vector<GroupStyle> &styles) {
    for (size_t i = 0; i < styles.size(); ++i) {
        if (styles[i].group == g) return (int)i;
    }
    return -1;
}

struct PlotPoint {
    double d_edge_deg;
    double xs_over_bh;
    double xs_over_bh_err;
    double bh_over_km15;
};

static void make_normalization_plots(const std::string &out_dir,
                                     const std::string &label,
                                     const std::string &helicity,
                                     const std::vector<PlotPoint> &pts) {
    ensure_output_dir_or_throw(out_dir);

    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);

    const std::string label_tag = sanitize_for_filename(label);
    const std::string hel_tag   = sanitize_for_filename(helicity);

    const double x_min = 0.1;
    double x_max = x_min * 10.0;
    for (size_t i = 0; i < pts.size(); ++i) {
        if (std::isfinite(pts[i].d_edge_deg) && pts[i].d_edge_deg > x_max) {
            x_max = pts[i].d_edge_deg;
        }
    }
    x_max *= 1.10;

    // 1D: xs/BH vs d_edge, grouped by BH/KM15 bins (legacy behavior)
    {
        const std::vector<GroupStyle> styles = make_group_styles();

        std::vector<TGraphErrors*> graphs(styles.size(), (TGraphErrors*)nullptr);
        std::vector<int> n_in(styles.size(), 0);

        for (size_t si = 0; si < styles.size(); ++si) {
            graphs[si] = new TGraphErrors();
            graphs[si]->SetName(Form("gr_%zu", si));
            graphs[si]->SetMarkerStyle(styles[si].marker_style);
            graphs[si]->SetMarkerColor(styles[si].marker_color);
            graphs[si]->SetLineColor(styles[si].marker_color);
            graphs[si]->SetLineWidth(1);
        }

        for (size_t i = 0; i < pts.size(); ++i) {
            const BhKmGroup g = categorize_bh_over_km15(pts[i].bh_over_km15);
            const int si = style_index_for_group(g, styles);
            if (si < 0) continue;

            const int n = n_in[(size_t)si];
            graphs[(size_t)si]->SetPoint(n, pts[i].d_edge_deg, pts[i].xs_over_bh);
            graphs[(size_t)si]->SetPointError(n, 0.0, pts[i].xs_over_bh_err);
            n_in[(size_t)si] += 1;
        }

        TCanvas *c = new TCanvas("c_norm_1d", "c_norm_1d", 950, 700);
        c->SetLogx(1);
        c->SetLeftMargin(0.14);
        c->SetRightMargin(0.05);
        c->SetBottomMargin(0.12);
        c->SetTopMargin(0.10);

        int first_nonempty = -1;
        for (size_t si = 0; si < graphs.size(); ++si) {
            if (graphs[si] && graphs[si]->GetN() > 0) {
                first_nonempty = (int)si;
                break;
            }
        }

        if (first_nonempty >= 0) {
            graphs[(size_t)first_nonempty]->Draw("AP");
            graphs[(size_t)first_nonempty]->GetXaxis()->SetTitle("d_edge (deg)");
            graphs[(size_t)first_nonempty]->GetYaxis()->SetTitle("xs/BH");
            graphs[(size_t)first_nonempty]->GetXaxis()->SetLimits(x_min, x_max);
            graphs[(size_t)first_nonempty]->GetYaxis()->SetRangeUser(0.0, 3.0);

            TLatex tl;
            tl.SetNDC(kTRUE);
            tl.SetTextSize(0.040);
            tl.DrawLatex(0.14, 0.93, Form("Normalization study (1D): %s, %s", label.c_str(), helicity.c_str()));

            TLine *l = new TLine(x_min, 1.0, x_max, 1.0);
            l->SetLineStyle(2);
            l->SetLineColor(kGray+2);
            l->SetLineWidth(2);
            l->Draw("same");

            for (size_t si = 0; si < graphs.size(); ++si) {
                if ((int)si == first_nonempty) continue;
                if (graphs[si] && graphs[si]->GetN() > 0) {
                    graphs[si]->Draw("P same");
                }
            }

            TLegend *leg = new TLegend(0.16, 0.66, 0.70, 0.88);
            leg->SetBorderSize(1);
            leg->SetFillStyle(1001);
            leg->SetFillColor(kWhite);
            leg->SetTextSize(0.028);

            for (size_t si = 0; si < graphs.size(); ++si) {
                if (graphs[si] && graphs[si]->GetN() > 0) {
                    leg->AddEntry(graphs[si], styles[si].label.c_str(), "p");
                }
            }
            leg->Draw();

            const std::string out_png =
                out_dir + "/norm_1d_xs_over_bh_vs_dedge_" + label_tag + "_" + hel_tag + ".png";
            c->SaveAs(out_png.c_str());
        } else {
            std::cerr << "[overall_norm] WARNING: no points available for 1D plot.\n";
        }

        for (size_t si = 0; si < graphs.size(); ++si) {
            delete graphs[si];
        }
        delete c;
    }

    // 2D: BH/KM15 vs d_edge, colored by xs/BH (legacy behavior)
    {
        TGraph2D *g2 = new TGraph2D();
        g2->SetName("g2_norm");

        int n = 0;
        for (size_t i = 0; i < pts.size(); ++i) {
            if (!std::isfinite(pts[i].d_edge_deg)) continue;
            if (!std::isfinite(pts[i].bh_over_km15)) continue;
            if (!std::isfinite(pts[i].xs_over_bh)) continue;

            g2->SetPoint(n, pts[i].d_edge_deg, pts[i].bh_over_km15, pts[i].xs_over_bh);
            n += 1;
        }

        if (g2->GetN() > 0) {
            TCanvas *c = new TCanvas("c_norm_2d", "c_norm_2d", 980, 720);
            c->SetLogx(1);
            c->SetLeftMargin(0.14);
            c->SetRightMargin(0.16);
            c->SetBottomMargin(0.12);
            c->SetTopMargin(0.10);

            g2->Draw("PCOLZ");

            if (g2->GetXaxis()) g2->GetXaxis()->SetLimits(x_min, x_max);

            g2->SetTitle("");
            if (g2->GetXaxis()) g2->GetXaxis()->SetTitle("d_edge (deg)");
            if (g2->GetYaxis()) g2->GetYaxis()->SetTitle("BH/KM15");
            if (g2->GetZaxis()) g2->GetZaxis()->SetTitle("xs/BH");

            TLatex tl;
            tl.SetNDC(kTRUE);
            tl.SetTextSize(0.040);
            tl.DrawLatex(0.14, 0.93, Form("Normalization study (2D): %s, %s", label.c_str(), helicity.c_str()));

            const std::string out_png =
                out_dir + "/norm_2d_bh_over_km15_vs_dedge_color_xs_over_bh_" + label_tag + "_" + hel_tag + ".png";
            c->SaveAs(out_png.c_str());

            delete c;
        } else {
            std::cerr << "[overall_norm] WARNING: no points available for 2D plot.\n";
        }

        delete g2;
    }
}

// -----------------------------------------------------------------------------
// New grid plotting: per (xB, Q2, |t|) bin fit xs(phi) with harmonic model,
// evaluate at phi0=0.001 deg, compare to BH(phi0), print & annotate.
// -----------------------------------------------------------------------------
struct PhiPointRow {
    double phi_deg;
    double xs;
    double xs_stat;
};

struct BinEdges {
    double min;
    double max;
};

static bool edges_equal(const BinEdges &a, const BinEdges &b) {
    // CSV edges are strings; we compare doubles here. Require exact equality
    // after atof-based parsing is not reliable; instead compare with a strict
    // tolerance small enough for typical bin edge formatting.
    const double tol = 1e-12;
    return (std::fabs(a.min - b.min) <= tol) && (std::fabs(a.max - b.max) <= tol);
}

static double eval_fit_uncertainty_at_phi0(const TFitResultPtr &res, const TF1 *f, double phi0_deg) {
    if (!res.Get()) return std::numeric_limits<double>::quiet_NaN();
    if (!f) return std::numeric_limits<double>::quiet_NaN();

    const TMatrixDSym cov = res->GetCovarianceMatrix();
    if (cov.GetNrows() < 4) return std::numeric_limits<double>::quiet_NaN();

    const double A = f->GetParameter(0);
    const double B = f->GetParameter(1);
    const double C = f->GetParameter(2);
    const double D = f->GetParameter(3);

    const double phi_rad = phi0_deg * TMath::DegToRad();
    const double c1 = std::cos(1.0 * phi_rad);
    const double c2 = std::cos(2.0 * phi_rad);
    const double c3 = std::cos(3.0 * phi_rad);

    // F = A*(1 + B*c1 + C*c2 + D*c3)
    const double dFdA = (1.0 + B * c1 + C * c2 + D * c3);
    const double dFdB = A * c1;
    const double dFdC = A * c2;
    const double dFdD = A * c3;

    const double g[4] = {dFdA, dFdB, dFdC, dFdD};

    double var = 0.0;
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            var += g[i] * cov(i, j) * g[j];
        }
    }

    if (!(var >= 0.0) || !std::isfinite(var)) return std::numeric_limits<double>::quiet_NaN();
    return std::sqrt(var);
}

static void draw_grid_for_xb_bin(const std::string &out_dir,
                                 const std::string &label,
                                 const std::string &helicity,
                                 double Ebeam,
                                 Helicity hel,
                                 const BinEdges &xb_edges,
                                 const std::vector<BinEdges> &q2_bins,
                                 const std::vector<BinEdges> &t_bins,
                                 const std::map<std::string, std::vector<PhiPointRow> > &phi_points_by_key,
                                 const std::map<std::string, double> &xbavg_by_key,
                                 const std::map<std::string, double> &q2avg_by_key,
                                 const std::map<std::string, double> &tavg_by_key) {
    ensure_output_dir_or_throw(out_dir);

    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);

    const int nrows = (int)q2_bins.size();
    const int ncols = (int)t_bins.size();
    if (nrows <= 0 || ncols <= 0) {
        return;
    }

    const int W = 420 * ncols + 200;
    const int H = 320 * nrows + 220;

    const std::string label_tag = sanitize_for_filename(label);
    const std::string hel_tag   = sanitize_for_filename(helicity);

    const std::string xb_tag =
        Form("xb_%0.4f_%0.4f", xb_edges.min, xb_edges.max);

    TCanvas *c = new TCanvas(Form("c_grid_%s_%s_%s", xb_tag.c_str(), label_tag.c_str(), hel_tag.c_str()),
                             "c_grid",
                             W, H);

    c->SetLeftMargin(0.08);
    c->SetRightMargin(0.02);
    c->SetBottomMargin(0.08);
    c->SetTopMargin(0.10);

    c->Divide(ncols, nrows, 0.001, 0.001);

    // Global title
    c->cd(1);
    TLatex title;
    title.SetNDC(kTRUE);
    title.SetTextSize(0.060);
    title.DrawLatex(0.10, 0.97,
        Form("Normalization grids: %s, %s   xB in [%0.3f, %0.3f]   Ebeam=%0.1f",
             label.c_str(), helicity.c_str(), xb_edges.min, xb_edges.max, Ebeam));

    // Iterate pads row-major: q2 row, t col
    for (int r = 0; r < nrows; ++r) {
        for (int ccol = 0; ccol < ncols; ++ccol) {
            const int ipad = r * ncols + ccol + 1;
            c->cd(ipad);

            gPad->SetLeftMargin(0.14);
            gPad->SetRightMargin(0.04);
            gPad->SetBottomMargin(0.14);
            gPad->SetTopMargin(0.18);

            const BinEdges q2e = q2_bins[(size_t)r];
            const BinEdges te  = t_bins[(size_t)ccol];

            // Build key using the same string formatting style as CSV fields is not safe here.
            // Instead, scan phi_points_by_key and match by xb/q2/t avg maps is not allowed (implicit).
            // So we require exact bin-edge matching by reconstructing key strings from numeric edges
            // using a deterministic format consistent with typical CSV.
            //
            // IMPORTANT: If your CSV uses different formatting (e.g. "0.1" vs "0.10"),
            // you must keep the CSV strings and use those. This code avoids fallbacks.
            //
            // Therefore: we do NOT synthesize keys here. We search for matching bins by comparing
            // parsed numeric edges stored in the key maps built in print_bh_normalization_study,
            // and we store a direct lookup map there. That is handled upstream by `phi_points_by_key`
            // containing only exact keys; here we search by numeric edges via a separate map.
            //
            // In this implementation, we pass only the `phi_points_by_key` keyed by the canonical
            // string key and additionally pass avg maps; for finding the correct key we will do a
            // deterministic scan over keys and compare the bin edges encoded in the key.
            //
            // This is explicit and deterministic (O(Nkeys)), and avoids silent formatting hacks.

            std::string matched_key;
            bool found = false;

            for (std::map<std::string, std::vector<PhiPointRow> >::const_iterator it = phi_points_by_key.begin();
                 it != phi_points_by_key.end(); ++it) {

                const std::string &key = it->first;

                // key format: xbmin|xbmax|q2min|q2max|tmin|tmax
                std::vector<std::string> parts;
                {
                    std::string tok;
                    for (size_t i = 0; i < key.size(); ++i) {
                        const char ch = key[i];
                        if (ch == '|') {
                            parts.push_back(tok);
                            tok.clear();
                        } else {
                            tok.push_back(ch);
                        }
                    }
                    parts.push_back(tok);
                }

                if (parts.size() != 6) continue;

                const double xbmin = std::atof(parts[0].c_str());
                const double xbmax = std::atof(parts[1].c_str());
                const double q2min = std::atof(parts[2].c_str());
                const double q2max = std::atof(parts[3].c_str());
                const double tmin  = std::atof(parts[4].c_str());
                const double tmax  = std::atof(parts[5].c_str());

                const BinEdges xb_test = {xbmin, xbmax};
                const BinEdges q2_test = {q2min, q2max};
                const BinEdges t_test  = {tmin,  tmax};

                if (!edges_equal(xb_test, xb_edges)) continue;
                if (!edges_equal(q2_test, q2e)) continue;
                if (!edges_equal(t_test, te)) continue;

                matched_key = key;
                found = true;
                break;
            }

            TLatex tl;
            tl.SetNDC(kTRUE);
            tl.SetTextSize(0.055);

            tl.DrawLatex(0.12, 0.92, Form("Q2 in [%0.3f,%0.3f]   |t| in [%0.3f,%0.3f]",
                                          q2e.min, q2e.max, te.min, te.max));

            if (!found) {
                tl.DrawLatex(0.12, 0.70, "No bin");
                continue;
            }

            const std::vector<PhiPointRow> &pts = phi_points_by_key.find(matched_key)->second;
            if (pts.empty()) {
                tl.DrawLatex(0.12, 0.70, "No points");
                continue;
            }

            // Build graph
            TGraphErrors *gr = new TGraphErrors();
            gr->SetName(Form("gr_%d_%d", r, ccol));
            gr->SetMarkerStyle(20);
            gr->SetMarkerSize(0.9);
            gr->SetLineWidth(1);

            int n = 0;
            for (size_t i = 0; i < pts.size(); ++i) {
                if (!std::isfinite(pts[i].phi_deg)) continue;
                if (!std::isfinite(pts[i].xs)) continue;
                if (!(pts[i].xs > 0.0)) continue;

                gr->SetPoint(n, pts[i].phi_deg, pts[i].xs);
                const double ey = (std::isfinite(pts[i].xs_stat) && pts[i].xs_stat >= 0.0) ? pts[i].xs_stat : 0.0;
                gr->SetPointError(n, 0.0, ey);
                n += 1;
            }

            if (gr->GetN() <= 0) {
                tl.DrawLatex(0.12, 0.70, "No valid points");
                delete gr;
                continue;
            }

            // Draw axes and points
            gr->Draw("AP");
            gr->GetXaxis()->SetTitle("#phi (deg)");
            gr->GetYaxis()->SetTitle("xs");
            gr->GetXaxis()->SetLimits(0.0, 360.0);

            // Set a reasonable y-range
            double ymax = 0.0;
            for (int i = 0; i < gr->GetN(); ++i) {
                double x, y;
                gr->GetPoint(i, x, y);
                if (std::isfinite(y) && y > ymax) ymax = y;
            }
            if (ymax <= 0.0) ymax = 1.0;
            gr->GetYaxis()->SetRangeUser(0.0, 1.25 * ymax);

            // Fit harmonic model if enough points
            // 4 parameters -> need at least 4 points, but realistically >= 5
            const int npts = gr->GetN();
            if (npts < 5) {
                tl.DrawLatex(0.12, 0.80, Form("Too few points to fit (N=%d)", npts));
                delete gr;
                continue;
            }

            TF1 *f = new TF1(Form("f_%d_%d", r, ccol), xs_harmonic_model_deg, 0.0, 360.0, 4);
            f->SetParNames("A", "B", "C", "D");

            // Initial values: A ~ mean y, others 0
            double sumy = 0.0;
            for (int i = 0; i < npts; ++i) {
                double x, y;
                gr->GetPoint(i, x, y);
                sumy += y;
            }
            const double A0 = (npts > 0) ? (sumy / (double)npts) : 1.0;

            f->SetParameter(0, (A0 > 0.0 ? A0 : 1.0));
            f->SetParameter(1, 0.0);
            f->SetParameter(2, 0.0);
            f->SetParameter(3, 0.0);

            // Quiet fit but keep result; draw the function
            const TFitResultPtr res = gr->Fit(f, "QS");

            const int fit_status = res.Get() ? res->Status() : 999;
            if (fit_status != 0) {
                tl.DrawLatex(0.12, 0.80, "Fit failed");
                delete f;
                delete gr;
                continue;
            }

            const double A  = f->GetParameter(0);
            const double dA = f->GetParError(0);

            // Evaluate xs_fit at phi0 and uncertainty (full covariance)
            const double xs_phi0 = f->Eval(kPhi0Deg);
            const double dxs_phi0 = eval_fit_uncertainty_at_phi0(res, f, kPhi0Deg);

            // Model BH and KM15 at phi0
            // t input for model functions is positive |t|
            double xb_c = std::numeric_limits<double>::quiet_NaN();
            double q2_c = std::numeric_limits<double>::quiet_NaN();
            double t_c  = std::numeric_limits<double>::quiet_NaN();

            {
                std::map<std::string, double>::const_iterator itx = xbavg_by_key.find(matched_key);
                std::map<std::string, double>::const_iterator itq = q2avg_by_key.find(matched_key);
                std::map<std::string, double>::const_iterator itt = tavg_by_key.find(matched_key);
                if (itx != xbavg_by_key.end()) xb_c = itx->second;
                if (itq != q2avg_by_key.end()) q2_c = itq->second;
                if (itt != tavg_by_key.end()) t_c  = itt->second;
            }

            double bh_phi0 = std::numeric_limits<double>::quiet_NaN();
            double km_phi0 = std::numeric_limits<double>::quiet_NaN();
            double bh_over_km_phi0 = std::numeric_limits<double>::quiet_NaN();
            double xs_over_bh_phi0 = std::numeric_limits<double>::quiet_NaN();

            if (finite_pos(xb_c) && finite_pos(q2_c) && finite_pos(t_c)) {
                bh_phi0 = vgg_bh_only(xb_c, q2_c, t_c, kPhi0Deg, Ebeam);
                km_phi0 = km15_xs(xb_c, q2_c, t_c, kPhi0Deg, Ebeam, hel);

                if (finite_pos(km_phi0) && finite_pos(bh_phi0)) {
                    bh_over_km_phi0 = bh_phi0 / km_phi0;
                }
                if (finite_pos(bh_phi0) && std::isfinite(xs_phi0)) {
                    xs_over_bh_phi0 = xs_phi0 / bh_phi0;
                }
            }

            // Annotate (no "Fit:" line)
            TLatex ann;
            ann.SetNDC(kTRUE);
            ann.SetTextSize(0.055);

            ann.DrawLatex(0.12, 0.82, Form("A = %0.6g +/- %0.3g", A, dA));

            if (std::isfinite(bh_phi0)) {
                ann.DrawLatex(0.12, 0.74, Form("BH(#phi=%.3g deg) = %0.6g", kPhi0Deg, bh_phi0));
            } else {
                ann.DrawLatex(0.12, 0.74, Form("BH(#phi=%.3g deg) = N/A", kPhi0Deg));
            }

            if (std::isfinite(km_phi0)) {
                ann.DrawLatex(0.12, 0.66, Form("KM15(#phi=%.3g deg) = %0.6g", kPhi0Deg, km_phi0));
            } else {
                ann.DrawLatex(0.12, 0.66, Form("KM15(#phi=%.3g deg) = N/A", kPhi0Deg));
            }

            if (std::isfinite(bh_over_km_phi0)) {
                ann.DrawLatex(0.12, 0.58, Form("BH/KM15 at #phi=%.3g deg: %0.6g", kPhi0Deg, bh_over_km_phi0));
            } else {
                ann.DrawLatex(0.12, 0.58, Form("BH/KM15 at #phi=%.3g deg: N/A", kPhi0Deg));
            }

            if (std::isfinite(xs_phi0)) {
                if (std::isfinite(dxs_phi0)) {
                    ann.DrawLatex(0.12, 0.50, Form("xs_fit(#phi=%.3g deg) = %0.6g +/- %0.3g", kPhi0Deg, xs_phi0, dxs_phi0));
                } else {
                    ann.DrawLatex(0.12, 0.50, Form("xs_fit(#phi=%.3g deg) = %0.6g", kPhi0Deg, xs_phi0));
                }
            } else {
                ann.DrawLatex(0.12, 0.50, Form("xs_fit(#phi=%.3g deg) = N/A", kPhi0Deg));
            }

            if (std::isfinite(xs_over_bh_phi0)) {
                ann.DrawLatex(0.12, 0.42, Form("xs_fit/BH at #phi=%.3g deg: %0.6g", kPhi0Deg, xs_over_bh_phi0));
            } else {
                ann.DrawLatex(0.12, 0.42, Form("xs_fit/BH at #phi=%.3g deg: N/A", kPhi0Deg));
            }

            delete f;
            delete gr;
        }
    }

    const std::string out_png =
        out_dir + "/norm_grid_" + xb_tag + "_" + label_tag + "_" + hel_tag + ".png";
    c->SaveAs(out_png.c_str());

    delete c;
}

} // end anonymous namespace

bool print_bh_normalization_study(const std::string &csv_path,
                                  const std::string &label,
                                  const std::string &helicity) {
    try {
        std::ifstream ifs(csv_path.c_str());
        if (!ifs) {
            std::cerr << "[overall_norm] ERROR: cannot open CSV: " << csv_path << "\n";
            return false;
        }

        std::vector<std::string> lines;
        std::string line;
        while (std::getline(ifs, line)) {
            lines.push_back(line);
        }
        ifs.close();

        if (lines.empty()) {
            std::cerr << "[overall_norm] ERROR: CSV is empty: " << csv_path << "\n";
            return false;
        }

        std::vector<std::string> header = split_csv_line(lines[0]);

        const int c_xbmin = require_col(header, "xBmin");
        const int c_xbmax = require_col(header, "xBmax");
        const int c_q2min = require_col(header, "Q2min");
        const int c_q2max = require_col(header, "Q2max");
        const int c_tmin  = require_col(header, "t_abs_min");
        const int c_tmax  = require_col(header, "t_abs_max");

        const std::string col_xbavg  = "xBavg, " + label;
        const std::string col_q2avg  = "Q2avg, " + label;
        const std::string col_tavg   = "t_abs_avg, " + label;
        const std::string col_phiavg = "phiavg, " + label;

        const int c_xbavg  = require_col(header, col_xbavg);
        const int c_q2avg  = require_col(header, col_q2avg);
        const int c_tavg   = require_col(header, col_tavg);
        const int c_phiavg = require_col(header, col_phiavg);

        const std::string col_xs =
            "cross sections, ep->epg, exp, " + label + ", " + helicity;
        const int c_xs = require_col(header, col_xs);

        const double Ebeam = beam_energy_for_label(label);
        const Helicity hel = helicity_from_string(helicity);

        std::cout << "============================================================\n";
        std::cout << "[overall_norm] BH normalization study\n";
        std::cout << "[overall_norm] CSV      : " << csv_path << "\n";
        std::cout << "[overall_norm] label    : " << label << "\n";
        std::cout << "[overall_norm] helicity : " << helicity << "\n";
        std::cout << "[overall_norm] Ebeam    : " << Ebeam << "\n";
        std::cout << "------------------------------------------------------------\n";

        // Legacy: pick best (closest-to-edge) row per kinematic bin
        std::map<std::string, BestPhiRow> best;

        // New: keep all phi points per kinematic bin for harmonic fits
        std::map<std::string, std::vector<PhiPointRow> > phi_points_by_key;

        // For grid diagnostics and BH/KM15 at phi0, keep per-bin central values
        std::map<std::string, double> xbavg_by_key;
        std::map<std::string, double> q2avg_by_key;
        std::map<std::string, double> tavg_by_key;

        // Track unique xB bins, and within each xB bin the unique Q2 and |t| bins
        std::vector<BinEdges> xb_bins;
        std::map<std::string, BinEdges> xb_edges_by_key;
        std::map<std::string, BinEdges> q2_edges_by_key;
        std::map<std::string, BinEdges> t_edges_by_key;

        size_t n_rows = (lines.size() > 0 ? lines.size() - 1 : 0);
        for (size_t r = 1; r < lines.size(); ++r) {
            if (lines[r].empty()) continue;

            std::vector<std::string> fields = split_csv_line(lines[r]);
            if (fields.size() != header.size()) continue;

            const std::string xbmin_s = trim(unquote(fields[c_xbmin]));
            const std::string xbmax_s = trim(unquote(fields[c_xbmax]));
            const std::string q2min_s = trim(unquote(fields[c_q2min]));
            const std::string q2max_s = trim(unquote(fields[c_q2max]));
            const std::string tmin_s  = trim(unquote(fields[c_tmin]));
            const std::string tmax_s  = trim(unquote(fields[c_tmax]));

            const std::string key =
                cell_key_for_kin_bin(xbmin_s, xbmax_s, q2min_s, q2max_s, tmin_s, tmax_s);

            const double xbmin = std::atof(xbmin_s.c_str());
            const double xbmax = std::atof(xbmax_s.c_str());
            const double q2min = std::atof(q2min_s.c_str());
            const double q2max = std::atof(q2max_s.c_str());
            const double tmin  = std::atof(tmin_s.c_str());
            const double tmax  = std::atof(tmax_s.c_str());

            xb_edges_by_key[key] = BinEdges{xbmin, xbmax};
            q2_edges_by_key[key] = BinEdges{q2min, q2max};
            t_edges_by_key[key]  = BinEdges{tmin,  tmax};

            // central values and phi
            const double xb_c  = std::atof(trim(unquote(fields[c_xbavg])).c_str());
            const double q2_c  = std::atof(trim(unquote(fields[c_q2avg])).c_str());
            const double t_c   = std::atof(trim(unquote(fields[c_tavg])).c_str());
            const double phi   = std::atof(trim(unquote(fields[c_phiavg])).c_str());

            if (!finite_pos(xb_c) || !finite_pos(q2_c) || !finite_pos(t_c) || !std::isfinite(phi)) {
                continue;
            }

            TripleCell xs = parse_tuple3(fields[c_xs]);
            if (!(xs.value > 0.0) || !std::isfinite(xs.value)) {
                continue;
            }

            xbavg_by_key[key] = xb_c;
            q2avg_by_key[key] = q2_c;
            tavg_by_key[key]  = t_c;

            // Save point for harmonic fit map
            {
                PhiPointRow pr;
                pr.phi_deg = phi;
                pr.xs      = xs.value;
                pr.xs_stat = xs.stat;
                phi_points_by_key[key].push_back(pr);
            }

            // Legacy "best edge point" selection
            const double dist = min_dist_to_0_or_360(phi);

            std::map<std::string, BestPhiRow>::iterator it = best.find(key);
            if (it == best.end() || dist < it->second.dist_to_edge) {
                BestPhiRow br;
                br.xb_c = xb_c;
                br.q2_c = q2_c;
                br.t_c  = t_c;
                br.phi_deg = phi;
                br.dist_to_edge = dist;
                br.xs = xs.value;
                br.xs_stat = xs.stat;
                br.csv_row_index = r;
                best[key] = br;
            }
        }

        if (best.empty()) {
            std::cerr << "[overall_norm] WARNING: no bins found with positive cross section values.\n";
            return true;
        }

        std::cout << "[overall_norm] Unique kinematic bins found: " << best.size()
                  << " (from " << n_rows << " CSV data rows)\n\n";

        // ---------------------------------------------------------------------
        // Legacy printout (closest-to-edge point) + legacy plots
        // ---------------------------------------------------------------------
        std::cout
            << std::setw(8)  << "xB"
            << std::setw(10) << "Q2"
            << std::setw(12) << "-t"
            << std::setw(10) << "d_edge"
            << std::setw(14) << "xs"
            << std::setw(14) << "xs/BH"
            << std::setw(14) << "BH/VGG"
            << std::setw(14) << "BH/KM15"
            << "\n";
        std::cout << std::string(8+10+12+10+14*4, '-') << "\n";

        std::vector<PlotPoint> plot_pts;
        plot_pts.reserve(best.size());

        // Weighted mean accumulator for xs/BH in BH/KM15 in [0.95,1.05]
        double sumw = 0.0;
        double sumwx = 0.0;
        int n_weighted_used = 0;
        int n_in_95_105_total = 0;

        for (std::map<std::string, BestPhiRow>::const_iterator it = best.begin();
             it != best.end(); ++it) {

            const BestPhiRow &br = it->second;

            const double xb   = br.xb_c;
            const double q2   = br.q2_c;
            const double tpos = br.t_c;
            const double phi  = br.phi_deg;

            const double bh  = vgg_bh_only(xb, q2, tpos, phi, Ebeam);
            const double vgg = vgg_xs(xb, q2, tpos, phi, Ebeam, hel);
            const double km  = km15_xs(xb, q2, tpos, phi, Ebeam, hel);

            double xs_over_bh  = 0.0;
            double xs_over_bh_err = 0.0;
            double bh_over_vgg = 0.0;
            double bh_over_km  = 0.0;

            if (finite_pos(bh)) {
                xs_over_bh = br.xs / bh;
                if (std::isfinite(br.xs_stat) && br.xs_stat >= 0.0) {
                    xs_over_bh_err = br.xs_stat / bh;
                }
            }
            if (finite_pos(vgg)) {
                bh_over_vgg = bh / vgg;
            }
            if (finite_pos(km)) {
                bh_over_km = bh / km;
            }

            const double tneg = -tpos;

            std::cout
                << std::setw(8)  << std::fixed << std::setprecision(3) << xb
                << std::setw(10) << std::fixed << std::setprecision(2) << q2
                << std::setw(12) << std::fixed << std::setprecision(3) << tneg
                << std::setw(10) << std::fixed << std::setprecision(1) << br.dist_to_edge
                << std::setw(14) << std::scientific << std::setprecision(3) << br.xs
                << std::setw(14) << std::fixed << std::setprecision(3) << xs_over_bh
                << std::setw(14) << std::fixed << std::setprecision(3) << bh_over_vgg
                << std::setw(14) << std::fixed << std::setprecision(3) << bh_over_km
                << "\n";

            if (std::isfinite(br.dist_to_edge) && br.dist_to_edge > 0.0 &&
                std::isfinite(xs_over_bh) && xs_over_bh >= 0.0 &&
                std::isfinite(bh_over_km) && bh_over_km > 0.0) {

                PlotPoint p;
                p.d_edge_deg      = br.dist_to_edge;
                p.xs_over_bh      = xs_over_bh;
                p.xs_over_bh_err  = xs_over_bh_err;
                p.bh_over_km15    = bh_over_km;
                plot_pts.push_back(p);
            }

            if (std::isfinite(bh_over_km) && bh_over_km >= 0.95 && bh_over_km <= 1.05) {
                n_in_95_105_total += 1;

                if (std::isfinite(xs_over_bh_err) && xs_over_bh_err > 0.0 &&
                    std::isfinite(xs_over_bh)) {

                    const double w = 1.0 / (xs_over_bh_err * xs_over_bh_err);
                    sumw  += w;
                    sumwx += w * xs_over_bh;
                    n_weighted_used += 1;
                }
            }
        }

        const std::string out_dir = "output/normalization_study";
        make_normalization_plots(out_dir, label, helicity, plot_pts);

        // Legacy weighted-mean summary
        std::cout << "\n";
        std::cout << "------------------------------------------------------------\n";
        std::cout << "[overall_norm] Weighted xs/BH for BH/KM15 in [0.95, 1.05]\n";
        std::cout << "[overall_norm] Points in range (total) : " << n_in_95_105_total << "\n";
        std::cout << "[overall_norm] Points used (err>0)     : " << n_weighted_used << "\n";

        if (sumw > 0.0) {
            const double mean = sumwx / sumw;
            const double err  = std::sqrt(1.0 / sumw);
            std::cout << "[overall_norm] Weighted mean xs/BH     : " << std::fixed << std::setprecision(6) << mean << "\n";
            std::cout << "[overall_norm] Weighted stat unc       : " << std::fixed << std::setprecision(6) << err  << "\n";
        } else {
            std::cout << "[overall_norm] Weighted mean xs/BH     : N/A (no usable uncertainties)\n";
        }
        std::cout << "------------------------------------------------------------\n";

        // ---------------------------------------------------------------------
        // New: per-bin harmonic fits evaluated at phi0 and printed (with BH/KM15)
        // ---------------------------------------------------------------------
        std::cout << "\n";
        std::cout << "------------------------------------------------------------\n";
        std::cout << "[overall_norm] Per-bin harmonic fits: xs(#phi) = A*(1 + B cos#phi + C cos2#phi + D cos3#phi)\n";
        std::cout << "[overall_norm] Evaluation at #phi = " << kPhi0Deg << " (deg)\n";
        std::cout
            << std::setw(8)  << "xB"
            << std::setw(10) << "Q2"
            << std::setw(12) << "-t"
            << std::setw(14) << "A"
            << std::setw(14) << "dA(stat)"
            << std::setw(14) << "BH(phi0)"
            << std::setw(14) << "KM15(phi0)"
            << std::setw(14) << "BH/KM15"
            << std::setw(14) << "xs_fit(phi0)"
            << std::setw(14) << "dxs_fit"
            << "\n";
        std::cout << std::string(8+10+12+14*7, '-') << "\n";

        // Build unique xB bins from keys
        {
            std::set<std::string> xb_seen;
            for (std::map<std::string, BinEdges>::const_iterator it = xb_edges_by_key.begin();
                 it != xb_edges_by_key.end(); ++it) {

                const std::string &key = it->first;
                std::map<std::string, BinEdges>::const_iterator itx = xb_edges_by_key.find(key);
                if (itx == xb_edges_by_key.end()) continue;

                const BinEdges xb = itx->second;
                const std::string tag = Form("%0.12g|%0.12g", xb.min, xb.max);

                if (xb_seen.find(tag) == xb_seen.end()) {
                    xb_seen.insert(tag);
                    xb_bins.push_back(xb);
                }
            }
            std::sort(xb_bins.begin(), xb_bins.end(),
                [](const BinEdges &a, const BinEdges &b) -> bool {
                    if (a.min < b.min) return true;
                    if (a.min > b.min) return false;
                    return a.max < b.max;
                });
        }

        // Create grids per xB bin
        const std::string grid_dir = out_dir + "/grids";
        ensure_output_dir_or_throw(grid_dir);

        for (size_t ix = 0; ix < xb_bins.size(); ++ix) {
            const BinEdges xbE = xb_bins[ix];

            // Collect Q2 and t bins for this xB
            std::vector<BinEdges> q2_bins;
            std::vector<BinEdges> t_bins;

            {
                std::set<std::string> q2_seen;
                std::set<std::string> t_seen;

                for (std::map<std::string, std::vector<PhiPointRow> >::const_iterator it = phi_points_by_key.begin();
                     it != phi_points_by_key.end(); ++it) {

                    const std::string &key = it->first;

                    std::map<std::string, BinEdges>::const_iterator itxb = xb_edges_by_key.find(key);
                    std::map<std::string, BinEdges>::const_iterator itq2 = q2_edges_by_key.find(key);
                    std::map<std::string, BinEdges>::const_iterator itt  = t_edges_by_key.find(key);

                    if (itxb == xb_edges_by_key.end() || itq2 == q2_edges_by_key.end() || itt == t_edges_by_key.end()) {
                        continue;
                    }

                    if (!edges_equal(itxb->second, xbE)) continue;

                    const BinEdges q2e = itq2->second;
                    const BinEdges te  = itt->second;

                    const std::string q2_tag = Form("%0.12g|%0.12g", q2e.min, q2e.max);
                    const std::string t_tag  = Form("%0.12g|%0.12g", te.min, te.max);

                    if (q2_seen.find(q2_tag) == q2_seen.end()) {
                        q2_seen.insert(q2_tag);
                        q2_bins.push_back(q2e);
                    }
                    if (t_seen.find(t_tag) == t_seen.end()) {
                        t_seen.insert(t_tag);
                        t_bins.push_back(te);
                    }
                }

                std::sort(q2_bins.begin(), q2_bins.end(),
                    [](const BinEdges &a, const BinEdges &b) -> bool {
                        if (a.min < b.min) return true;
                        if (a.min > b.min) return false;
                        return a.max < b.max;
                    });
                std::sort(t_bins.begin(), t_bins.end(),
                    [](const BinEdges &a, const BinEdges &b) -> bool {
                        if (a.min < b.min) return true;
                        if (a.min > b.min) return false;
                        return a.max < b.max;
                    });
            }

            // Print per-bin fit summary lines for bins in this xB (stdout)
            for (std::map<std::string, std::vector<PhiPointRow> >::const_iterator it = phi_points_by_key.begin();
                 it != phi_points_by_key.end(); ++it) {

                const std::string &key = it->first;

                if (xb_edges_by_key.find(key) == xb_edges_by_key.end()) continue;
                if (!edges_equal(xb_edges_by_key.find(key)->second, xbE)) continue;

                // Build a temporary graph and fit (stdout only)
                const std::vector<PhiPointRow> &pts = it->second;
                if (pts.size() < 5) continue;

                TGraphErrors gr;
                int n = 0;
                for (size_t i = 0; i < pts.size(); ++i) {
                    if (!std::isfinite(pts[i].phi_deg)) continue;
                    if (!std::isfinite(pts[i].xs) || !(pts[i].xs > 0.0)) continue;
                    gr.SetPoint(n, pts[i].phi_deg, pts[i].xs);
                    const double ey = (std::isfinite(pts[i].xs_stat) && pts[i].xs_stat >= 0.0) ? pts[i].xs_stat : 0.0;
                    gr.SetPointError(n, 0.0, ey);
                    n += 1;
                }
                if (gr.GetN() < 5) continue;

                TF1 f("f_tmp", xs_harmonic_model_deg, 0.0, 360.0, 4);
                f.SetParNames("A", "B", "C", "D");

                double sumy = 0.0;
                for (int i = 0; i < gr.GetN(); ++i) {
                    double x, y;
                    gr.GetPoint(i, x, y);
                    sumy += y;
                }
                const double A0 = sumy / (double)gr.GetN();
                f.SetParameter(0, (A0 > 0.0 ? A0 : 1.0));
                f.SetParameter(1, 0.0);
                f.SetParameter(2, 0.0);
                f.SetParameter(3, 0.0);

                const TFitResultPtr res = gr.Fit(&f, "QSN");
                if (!res.Get() || res->Status() != 0) continue;

                const double A  = f.GetParameter(0);
                const double dA = f.GetParError(0);

                const double xs_phi0  = f.Eval(kPhi0Deg);
                const double dxs_phi0 = eval_fit_uncertainty_at_phi0(res, &f, kPhi0Deg);

                const double xb_c = xbavg_by_key.find(key)->second;
                const double q2_c = q2avg_by_key.find(key)->second;
                const double t_c  = tavg_by_key.find(key)->second;

                const double bh_phi0 = vgg_bh_only(xb_c, q2_c, t_c, kPhi0Deg, Ebeam);
                const double km_phi0 = km15_xs(xb_c, q2_c, t_c, kPhi0Deg, Ebeam, hel);
                const double bh_over_km = (finite_pos(km_phi0) && finite_pos(bh_phi0)) ? (bh_phi0 / km_phi0) : std::numeric_limits<double>::quiet_NaN();

                const double tneg = -t_c;

                std::cout
                    << std::setw(8)  << std::fixed << std::setprecision(3) << xb_c
                    << std::setw(10) << std::fixed << std::setprecision(2) << q2_c
                    << std::setw(12) << std::fixed << std::setprecision(3) << tneg
                    << std::setw(14) << std::scientific << std::setprecision(6) << A
                    << std::setw(14) << std::scientific << std::setprecision(3) << dA
                    << std::setw(14) << std::scientific << std::setprecision(6) << bh_phi0
                    << std::setw(14) << std::scientific << std::setprecision(6) << km_phi0
                    << std::setw(14) << std::scientific << std::setprecision(6) << bh_over_km
                    << std::setw(14) << std::scientific << std::setprecision(6) << xs_phi0
                    << std::setw(14) << std::scientific << std::setprecision(3) << dxs_phi0
                    << "\n";
            }

            // Draw and save the grid canvas for this xB
            draw_grid_for_xb_bin(grid_dir, label, helicity, Ebeam, hel, xbE,
                                 q2_bins, t_bins, phi_points_by_key,
                                 xbavg_by_key, q2avg_by_key, tavg_by_key);
        }

        std::cout << "------------------------------------------------------------\n";
        std::cout << "[overall_norm] Done.\n";
        std::cout << "============================================================\n";
        return true;

    } catch (const std::exception &e) {
        std::cerr << "[overall_norm] FATAL: " << e.what() << "\n";
        return false;
    }
}