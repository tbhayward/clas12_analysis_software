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
#include <TPad.h>

#include <filesystem>

namespace {

static const double kMinDEdgeDeg = 0.001; // requested: 0.001 degrees

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
    int idx = find_col_optional(header, target);
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
        char first = s.front();
        char last  = s.back();
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
    // bin edges (also preserved as strings for deterministic keys / filenames)
    std::string xbmin_s, xbmax_s;
    std::string q2min_s, q2max_s;
    std::string tmin_s,  tmax_s;

    double xbmin, xbmax;
    double q2min, q2max;
    double tmin,  tmax;

    // centers + chosen phi row
    double xb_c;
    double q2_c;
    double t_c;   // positive |t|
    double phi_deg;
    double dist_to_edge;

    // measured cross section (exp)
    double xs;
    double xs_stat;
    size_t csv_row_index;

    // derived ratios (computed after selection)
    double xs_over_bh;
    double xs_over_bh_err;
    double bh_over_vgg;
    double bh_over_km15;

    BestPhiRow() :
        xbmin_s(""), xbmax_s(""),
        q2min_s(""), q2max_s(""),
        tmin_s(""),  tmax_s(""),
        xbmin(std::numeric_limits<double>::quiet_NaN()),
        xbmax(std::numeric_limits<double>::quiet_NaN()),
        q2min(std::numeric_limits<double>::quiet_NaN()),
        q2max(std::numeric_limits<double>::quiet_NaN()),
        tmin(std::numeric_limits<double>::quiet_NaN()),
        tmax(std::numeric_limits<double>::quiet_NaN()),
        xb_c(std::numeric_limits<double>::quiet_NaN()),
        q2_c(std::numeric_limits<double>::quiet_NaN()),
        t_c(std::numeric_limits<double>::quiet_NaN()),
        phi_deg(std::numeric_limits<double>::quiet_NaN()),
        dist_to_edge(std::numeric_limits<double>::infinity()),
        xs(0.0),
        xs_stat(0.0),
        csv_row_index(0),
        xs_over_bh(0.0),
        xs_over_bh_err(0.0),
        bh_over_vgg(0.0),
        bh_over_km15(0.0) {}
};

static double min_dist_to_0_or_360(double phi_deg) {
    double d0   = std::fabs(phi_deg - 0.0);
    double d360 = std::fabs(phi_deg - 360.0);
    double d = std::min(d0, d360);
    if (!std::isfinite(d)) return std::numeric_limits<double>::infinity();
    if (d < kMinDEdgeDeg) d = kMinDEdgeDeg; // explicit floor to support log-x plots
    return d;
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

// 5% bands around 1.00, with symmetric paired bins.
enum class BhKmGroup {
    G_95_105 = 0,
    G_90_95_OR_105_110,
    G_85_90_OR_110_115,
    G_OUTSIDE_15PCT,
    G_INVALID
};

static BhKmGroup categorize_bh_over_km15(double r) {
    if (!std::isfinite(r) || r <= 0.0) return BhKmGroup::G_INVALID;

    // Central 5% band: inclusive
    if (r >= 0.95 && r <= 1.05) return BhKmGroup::G_95_105;

    // Next band: [0.90, 0.95) and (1.05, 1.10]
    if ((r >= 0.90 && r < 0.95) || (r > 1.05 && r <= 1.10)) return BhKmGroup::G_90_95_OR_105_110;

    // Next band: (0.85, 0.90) and (1.10, 1.15)
    // Endpoints 0.85 and 1.15 are treated as "outside 15%" (inclusive outside band).
    if ((r > 0.85 && r < 0.90) || (r > 1.10 && r < 1.15)) return BhKmGroup::G_85_90_OR_110_115;

    // Outside 15% inclusive
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

    const double x_min = kMinDEdgeDeg;
    double x_max = x_min * 10.0;
    for (size_t i = 0; i < pts.size(); ++i) {
        if (std::isfinite(pts[i].d_edge_deg) && pts[i].d_edge_deg > x_max) {
            x_max = pts[i].d_edge_deg;
        }
    }
    x_max *= 1.10;

    // -------------------------------------------------------------------------
    // 1D: xs/BH vs d_edge, grouped by BH/KM15 bins
    // -------------------------------------------------------------------------
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

            const std::string out_png = out_dir + "/norm_1d_xs_over_bh_vs_dedge_" + label_tag + "_" + hel_tag + ".png";
            c->SaveAs(out_png.c_str());

            delete leg;
            delete l;
        } else {
            std::cerr << "[overall_norm] WARNING: no points available for 1D plot.\n";
        }

        for (size_t si = 0; si < graphs.size(); ++si) {
            delete graphs[si];
        }
        delete c;
    }

    // -------------------------------------------------------------------------
    // 2D: BH/KM15 vs d_edge, colored by xs/BH
    // -------------------------------------------------------------------------
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

            const std::string out_png = out_dir + "/norm_2d_bh_over_km15_vs_dedge_color_xs_over_bh_" + label_tag + "_" + hel_tag + ".png";
            c->SaveAs(out_png.c_str());

            delete c;
        } else {
            std::cerr << "[overall_norm] WARNING: no points available for 2D plot.\n";
        }

        delete g2;
    }
}

struct BinDef {
    std::string key;
    double minv;
    double maxv;
};

static bool bin_less_by_min(const BinDef &a, const BinDef &b) {
    if (a.minv < b.minv) return true;
    if (a.minv > b.minv) return false;
    return a.maxv < b.maxv;
}

static void make_xb_grid_canvases(const std::string &out_dir,
                                  const std::string &label,
                                  const std::string &helicity,
                                  const std::map<std::string, BestPhiRow> &best_rows) {
    ensure_output_dir_or_throw(out_dir);

    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);

    const std::vector<GroupStyle> styles = make_group_styles();

    // Group by xB bin (use the exact CSV edge strings as the canonical key).
    struct XbGroup {
        std::string xb_key;
        std::string xbmin_s;
        std::string xbmax_s;
        double xbmin;
        double xbmax;
        std::map<std::string, BinDef> q2bins_by_key;
        std::map<std::string, BinDef> tbins_by_key;
        std::map<std::string, std::map<std::string, const BestPhiRow*> > cell; // cell[q2key][tkey] -> row
    };

    std::map<std::string, XbGroup> groups;

    for (std::map<std::string, BestPhiRow>::const_iterator it = best_rows.begin();
         it != best_rows.end(); ++it) {

        const BestPhiRow &br = it->second;

        if (!(finite_pos(br.xbmin) && finite_pos(br.xbmax) && finite_pos(br.q2min) && finite_pos(br.q2max) &&
              finite_pos(br.tmin)  && finite_pos(br.tmax))) {
            continue;
        }

        const std::string xb_key = br.xbmin_s + "|" + br.xbmax_s;
        if (groups.find(xb_key) == groups.end()) {
            XbGroup g;
            g.xb_key  = xb_key;
            g.xbmin_s = br.xbmin_s;
            g.xbmax_s = br.xbmax_s;
            g.xbmin   = br.xbmin;
            g.xbmax   = br.xbmax;
            groups[xb_key] = g;
        }

        XbGroup &G = groups[xb_key];

        const std::string q2_key = br.q2min_s + "|" + br.q2max_s;
        if (G.q2bins_by_key.find(q2_key) == G.q2bins_by_key.end()) {
            BinDef b;
            b.key  = q2_key;
            b.minv = br.q2min;
            b.maxv = br.q2max;
            G.q2bins_by_key[q2_key] = b;
        }

        const std::string t_key = br.tmin_s + "|" + br.tmax_s;
        if (G.tbins_by_key.find(t_key) == G.tbins_by_key.end()) {
            BinDef b;
            b.key  = t_key;
            b.minv = br.tmin;
            b.maxv = br.tmax;
            G.tbins_by_key[t_key] = b;
        }

        G.cell[q2_key][t_key] = &br;
    }

    if (groups.empty()) {
        std::cerr << "[overall_norm] WARNING: no xB groups available for grid canvases.\n";
        return;
    }

    const std::string label_tag = sanitize_for_filename(label);
    const std::string hel_tag   = sanitize_for_filename(helicity);

    // For deterministic output order: sort by xbmin numeric, then xbmax numeric.
    std::vector<XbGroup> ordered;
    ordered.reserve(groups.size());
    for (std::map<std::string, XbGroup>::const_iterator it = groups.begin(); it != groups.end(); ++it) {
        ordered.push_back(it->second);
    }
    std::sort(ordered.begin(), ordered.end(),
              [](const XbGroup &a, const XbGroup &b) {
                  if (a.xbmin < b.xbmin) return true;
                  if (a.xbmin > b.xbmin) return false;
                  return a.xbmax < b.xbmax;
              });

    for (size_t gi = 0; gi < ordered.size(); ++gi) {
        const XbGroup &G = ordered[gi];

        // Build sorted q2 and t bin lists.
        std::vector<BinDef> q2bins;
        std::vector<BinDef> tbins;

        for (std::map<std::string, BinDef>::const_iterator it = G.q2bins_by_key.begin();
             it != G.q2bins_by_key.end(); ++it) {
            q2bins.push_back(it->second);
        }
        for (std::map<std::string, BinDef>::const_iterator it = G.tbins_by_key.begin();
             it != G.tbins_by_key.end(); ++it) {
            tbins.push_back(it->second);
        }

        std::sort(q2bins.begin(), q2bins.end(), bin_less_by_min);
        std::sort(tbins.begin(), tbins.end(), bin_less_by_min);

        const int nrows = (int)q2bins.size();
        const int ncols = (int)tbins.size();

        if (nrows <= 0 || ncols <= 0) continue;

        // Determine x_max for this canvas from the points present.
        const double x_min = kMinDEdgeDeg;
        double x_max = x_min * 10.0;

        double y_max = 3.0; // fixed for now (and consistent across subpads)
        (void)y_max;

        for (int ir = 0; ir < nrows; ++ir) {
            for (int ic = 0; ic < ncols; ++ic) {
                const std::string &q2_key = q2bins[(size_t)ir].key;
                const std::string &t_key  = tbins[(size_t)ic].key;

                std::map<std::string, std::map<std::string, const BestPhiRow*> >::const_iterator itq = G.cell.find(q2_key);
                if (itq == G.cell.end()) continue;
                std::map<std::string, const BestPhiRow*>::const_iterator itt = itq->second.find(t_key);
                if (itt == itq->second.end()) continue;

                const BestPhiRow *br = itt->second;
                if (!br) continue;
                if (std::isfinite(br->dist_to_edge) && br->dist_to_edge > x_max) {
                    x_max = br->dist_to_edge;
                }
            }
        }
        x_max *= 1.10;

        // Per-xB weighted mean (only within BH/KM15 in [0.95, 1.05], using xs/BH errors)
        double sumw = 0.0;
        double sumwx = 0.0;
        int n_in_95_105 = 0;
        int n_used = 0;

        for (int ir = 0; ir < nrows; ++ir) {
            for (int ic = 0; ic < ncols; ++ic) {
                const std::string &q2_key = q2bins[(size_t)ir].key;
                const std::string &t_key  = tbins[(size_t)ic].key;

                std::map<std::string, std::map<std::string, const BestPhiRow*> >::const_iterator itq = G.cell.find(q2_key);
                if (itq == G.cell.end()) continue;
                std::map<std::string, const BestPhiRow*>::const_iterator itt = itq->second.find(t_key);
                if (itt == itq->second.end()) continue;

                const BestPhiRow *br = itt->second;
                if (!br) continue;

                if (std::isfinite(br->bh_over_km15) && br->bh_over_km15 >= 0.95 && br->bh_over_km15 <= 1.05) {
                    n_in_95_105 += 1;
                    if (std::isfinite(br->xs_over_bh_err) && br->xs_over_bh_err > 0.0 &&
                        std::isfinite(br->xs_over_bh)) {
                        const double w = 1.0 / (br->xs_over_bh_err * br->xs_over_bh_err);
                        sumw  += w;
                        sumwx += w * br->xs_over_bh;
                        n_used += 1;
                    }
                }
            }
        }

        bool have_weighted = (sumw > 0.0);
        double wmean = 0.0;
        double werr  = 0.0;
        if (have_weighted) {
            wmean = sumwx / sumw;
            werr  = std::sqrt(1.0 / sumw);
        }

        const int W = 300 * ncols + 220;
        const int H = 260 * nrows + 240;

        const std::string cname = "c_norm_grid_xb_" + sanitize_for_filename(G.xb_key) + "_" + label_tag + "_" + hel_tag;
        TCanvas *c = new TCanvas(cname.c_str(), cname.c_str(), W, H);

        TPad *pTop  = new TPad("pTop",  "pTop",  0.0, 0.88, 1.0, 1.0);
        TPad *pGrid = new TPad("pGrid", "pGrid", 0.0, 0.00, 1.0, 0.88);

        pTop->SetMargin(0.0, 0.0, 0.0, 0.0);
        pGrid->SetMargin(0.0, 0.0, 0.0, 0.0);

        pTop->Draw();
        pGrid->Draw();

        pTop->cd();

        TLatex tl;
        tl.SetNDC(kTRUE);
        tl.SetTextSize(0.42);
        tl.DrawLatex(0.02, 0.62, Form("Normalization grid: %s, %s", label.c_str(), helicity.c_str()));
        tl.SetTextSize(0.40);
        tl.DrawLatex(0.02, 0.18, Form("xB bin: [%s, %s]", G.xbmin_s.c_str(), G.xbmax_s.c_str()));

        // Legend for BH/KM15 group colors in the top pad (dummy graphs).
        {
            TLegend *leg = new TLegend(0.52, 0.12, 0.98, 0.92);
            leg->SetBorderSize(1);
            leg->SetFillStyle(1001);
            leg->SetFillColor(kWhite);
            leg->SetTextSize(0.24);

            std::vector<TGraphErrors*> dummy;
            dummy.reserve(styles.size());
            for (size_t si = 0; si < styles.size(); ++si) {
                TGraphErrors *g = new TGraphErrors();
                g->SetMarkerStyle(styles[si].marker_style);
                g->SetMarkerColor(styles[si].marker_color);
                g->SetLineColor(styles[si].marker_color);
                g->SetLineWidth(1);
                dummy.push_back(g);
                leg->AddEntry(g, styles[si].label.c_str(), "p");
            }

            leg->Draw();

            // Weighted mean annotation (per-xB)
            if (have_weighted) {
                TLatex t2;
                t2.SetNDC(kTRUE);
                t2.SetTextSize(0.22);
                t2.DrawLatex(0.52, 0.02, Form("Weighted xs/BH (BH/KM15 in [0.95,1.05]): %.6f +/- %.6f (N=%d, used=%d)",
                                              wmean, werr, n_in_95_105, n_used));
            } else {
                TLatex t2;
                t2.SetNDC(kTRUE);
                t2.SetTextSize(0.22);
                t2.DrawLatex(0.52, 0.02, Form("Weighted xs/BH (BH/KM15 in [0.95,1.05]): N/A (N=%d, used=%d)",
                                              n_in_95_105, n_used));
            }

            // cleanup dummy graphs + legend
            for (size_t si = 0; si < dummy.size(); ++si) delete dummy[si];
            delete leg;
        }

        pGrid->cd();
        pGrid->Divide(ncols, nrows, 0.0, 0.0);

        for (int ir = 0; ir < nrows; ++ir) {
            for (int ic = 0; ic < ncols; ++ic) {
                const int ipad = ir * ncols + ic + 1;
                TPad *pad = (TPad*)pGrid->cd(ipad);
                if (!pad) continue;

                pad->SetLogx(1);
                pad->SetLeftMargin(0.16);
                pad->SetRightMargin(0.06);
                pad->SetBottomMargin(0.18);
                pad->SetTopMargin(0.16);
                pad->SetTicks(1, 1);

                // Frame
                TH1F *frame = (TH1F*)pad->DrawFrame(x_min, 0.0, x_max, 3.0);
                frame->SetTitle("");

                frame->GetXaxis()->SetTitle("d_edge (deg)");
                frame->GetYaxis()->SetTitle("xs/BH");

                frame->GetXaxis()->CenterTitle(true);
                frame->GetYaxis()->CenterTitle(true);

                frame->GetXaxis()->SetTitleSize(0.070);
                frame->GetYaxis()->SetTitleSize(0.070);

                frame->GetXaxis()->SetLabelSize(0.060);
                frame->GetYaxis()->SetLabelSize(0.060);

                frame->GetXaxis()->SetTitleOffset(1.10);
                frame->GetYaxis()->SetTitleOffset(1.15);

                // Reference line at 1.0
                TLine *l = new TLine(x_min, 1.0, x_max, 1.0);
                l->SetLineStyle(2);
                l->SetLineColor(kGray+2);
                l->SetLineWidth(2);
                l->Draw("same");

                const std::string &q2_key = q2bins[(size_t)ir].key;
                const std::string &t_key  = tbins[(size_t)ic].key;

                const BestPhiRow *br = nullptr;

                std::map<std::string, std::map<std::string, const BestPhiRow*> >::const_iterator itq = G.cell.find(q2_key);
                if (itq != G.cell.end()) {
                    std::map<std::string, const BestPhiRow*>::const_iterator itt = itq->second.find(t_key);
                    if (itt != itq->second.end()) br = itt->second;
                }

                // Bin label annotation
                TLatex tbin;
                tbin.SetNDC(kTRUE);
                tbin.SetTextSize(0.060);
                tbin.DrawLatex(0.12, 0.86, Form("Q2: [%s, %s] (GeV^2)", q2bins[(size_t)ir].key.substr(0, q2bins[(size_t)ir].key.find("|")).c_str(),
                                                q2bins[(size_t)ir].key.substr(q2bins[(size_t)ir].key.find("|")+1).c_str()));
                tbin.DrawLatex(0.12, 0.78, Form("|t|: [%s, %s] (GeV^2)", tbins[(size_t)ic].key.substr(0, tbins[(size_t)ic].key.find("|")).c_str(),
                                                tbins[(size_t)ic].key.substr(tbins[(size_t)ic].key.find("|")+1).c_str()));

                if (br && std::isfinite(br->xs_over_bh) && std::isfinite(br->dist_to_edge)) {
                    const BhKmGroup g = categorize_bh_over_km15(br->bh_over_km15);
                    const int si = style_index_for_group(g, styles);

                    int mstyle = 20;
                    int mcolor = kBlack;
                    if (si >= 0) {
                        mstyle = styles[(size_t)si].marker_style;
                        mcolor = styles[(size_t)si].marker_color;
                    }

                    TGraphErrors *gr = new TGraphErrors();
                    gr->SetPoint(0, br->dist_to_edge, br->xs_over_bh);
                    gr->SetPointError(0, 0.0, br->xs_over_bh_err);

                    gr->SetMarkerStyle(mstyle);
                    gr->SetMarkerColor(mcolor);
                    gr->SetLineColor(mcolor);
                    gr->SetLineWidth(1);

                    gr->Draw("P same");

                    TLatex tval;
                    tval.SetNDC(kTRUE);
                    tval.SetTextSize(0.060);
                    tval.DrawLatex(0.12, 0.68, Form("xs/BH = %.3f +/- %.3f", br->xs_over_bh, br->xs_over_bh_err));
                    tval.DrawLatex(0.12, 0.60, Form("BH/KM15 = %.3f", br->bh_over_km15));
                    tval.DrawLatex(0.12, 0.52, Form("d_edge = %.3f (deg)", br->dist_to_edge));

                    delete gr;
                } else {
                    TLatex tempty;
                    tempty.SetNDC(kTRUE);
                    tempty.SetTextSize(0.070);
                    tempty.DrawLatex(0.30, 0.40, "No point");
                }

                delete l;
            }
        }

        const std::string xbmin_tag = sanitize_for_filename(G.xbmin_s);
        const std::string xbmax_tag = sanitize_for_filename(G.xbmax_s);

        const std::string out_png =
            out_dir + "/norm_grid_xb_" + xbmin_tag + "_" + xbmax_tag + "_" + label_tag + "_" + hel_tag + ".png";

        c->SaveAs(out_png.c_str());

        delete pTop;
        delete pGrid;
        delete c;
    }
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

        std::map<std::string, BestPhiRow> best;

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

            const std::string key = cell_key_for_kin_bin(xbmin_s, xbmax_s, q2min_s, q2max_s, tmin_s, tmax_s);

            const double xbmin = std::atof(xbmin_s.c_str());
            const double xbmax = std::atof(xbmax_s.c_str());
            const double q2min = std::atof(q2min_s.c_str());
            const double q2max = std::atof(q2max_s.c_str());
            const double tmin  = std::atof(tmin_s.c_str());
            const double tmax  = std::atof(tmax_s.c_str());

            const double xb_c   = std::atof(trim(unquote(fields[c_xbavg])).c_str());
            const double q2_c   = std::atof(trim(unquote(fields[c_q2avg])).c_str());
            const double t_c    = std::atof(trim(unquote(fields[c_tavg])).c_str());
            const double phi    = std::atof(trim(unquote(fields[c_phiavg])).c_str());

            if (!finite_pos(xb_c) || !finite_pos(q2_c) || !finite_pos(t_c) || !std::isfinite(phi)) {
                continue;
            }

            TripleCell xs = parse_tuple3(fields[c_xs]);
            if (!(xs.value > 0.0) || !std::isfinite(xs.value)) {
                continue;
            }

            const double dist = min_dist_to_0_or_360(phi);

            std::map<std::string, BestPhiRow>::iterator it = best.find(key);
            if (it == best.end() || dist < it->second.dist_to_edge) {
                BestPhiRow br;

                br.xbmin_s = xbmin_s; br.xbmax_s = xbmax_s;
                br.q2min_s = q2min_s; br.q2max_s = q2max_s;
                br.tmin_s  = tmin_s;  br.tmax_s  = tmax_s;

                br.xbmin = xbmin; br.xbmax = xbmax;
                br.q2min = q2min; br.q2max = q2max;
                br.tmin  = tmin;  br.tmax  = tmax;

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

        // Weighted mean accumulator for BH/KM15 in [0.95, 1.05]
        double sumw = 0.0;
        double sumwx = 0.0;
        int n_weighted_used = 0;
        int n_in_95_105_total = 0;

        for (std::map<std::string, BestPhiRow>::iterator it = best.begin();
             it != best.end(); ++it) {

            BestPhiRow &br = it->second;

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

            br.xs_over_bh     = xs_over_bh;
            br.xs_over_bh_err = xs_over_bh_err;
            br.bh_over_vgg    = bh_over_vgg;
            br.bh_over_km15   = bh_over_km;

            const double tneg = -tpos;

            std::cout
                << std::setw(8)  << std::fixed << std::setprecision(3) << xb
                << std::setw(10) << std::fixed << std::setprecision(2) << q2
                << std::setw(12) << std::fixed << std::setprecision(3) << tneg
                << std::setw(10) << std::fixed << std::setprecision(3) << br.dist_to_edge
                << std::setw(14) << std::scientific << std::setprecision(3) << br.xs
                << std::setw(14) << std::fixed << std::setprecision(3) << xs_over_bh
                << std::setw(14) << std::fixed << std::setprecision(3) << bh_over_vgg
                << std::setw(14) << std::fixed << std::setprecision(3) << bh_over_km
                << "\n";

            if (std::isfinite(br.dist_to_edge) && br.dist_to_edge >= kMinDEdgeDeg &&
                std::isfinite(xs_over_bh) && xs_over_bh >= 0.0 &&
                std::isfinite(bh_over_km) && bh_over_km > 0.0) {

                PlotPoint p;
                p.d_edge_deg      = br.dist_to_edge;
                p.xs_over_bh      = xs_over_bh;
                p.xs_over_bh_err  = xs_over_bh_err;
                p.bh_over_km15    = bh_over_km;
                plot_pts.push_back(p);
            }

            // Weighted mean for BH/KM15 in [0.95, 1.05] inclusive
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

        // New: per-xB canvases with Q2 (rows) and |t| (cols)
        const std::string out_dir_grids = "output/normalization_study/grids";
        make_xb_grid_canvases(out_dir_grids, label, helicity, best);

        // Print weighted mean summary at the very end
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

        std::cout << "\n";
        std::cout << "[overall_norm] Done.\n";
        std::cout << "============================================================\n";
        return true;

    } catch (const std::exception &e) {
        std::cerr << "[overall_norm] FATAL: " << e.what() << "\n";
        return false;
    }
}