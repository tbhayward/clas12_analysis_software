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
#include <TGraph.h>
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

static std::string xb_key_for_edges(const std::string &xbmin_s,
                                    const std::string &xbmax_s) {
    return xbmin_s + "|" + xbmax_s;
}

static std::string q2_key_for_edges(const std::string &q2min_s,
                                    const std::string &q2max_s) {
    return q2min_s + "|" + q2max_s;
}

static std::string t_key_for_edges(const std::string &tmin_s,
                                   const std::string &tmax_s) {
    return tmin_s + "|" + tmax_s;
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

// Legacy grouping: 5% bands around 1.00, with symmetric paired bins.
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

struct PhiPoint {
    double phi_deg;
    double xs;
    double xs_stat;
    double dist_to_edge;

    PhiPoint() :
        phi_deg(std::numeric_limits<double>::quiet_NaN()),
        xs(0.0),
        xs_stat(0.0),
        dist_to_edge(std::numeric_limits<double>::infinity()) {}
};

struct KinBinData {
    // Edge strings (from CSV) used for deterministic bin identity and labeling
    std::string xbmin_s;
    std::string xbmax_s;
    std::string q2min_s;
    std::string q2max_s;
    std::string tmin_s;
    std::string tmax_s;

    // Numeric edges for sorting and pad labeling
    double xbmin;
    double xbmax;
    double q2min;
    double q2max;
    double tmin;
    double tmax;

    // Centers (from CSV)
    double xb_c;
    double q2_c;
    double t_c; // positive |t|

    std::vector<PhiPoint> points;

    KinBinData() :
        xbmin(0.0), xbmax(0.0),
        q2min(0.0), q2max(0.0),
        tmin(0.0),  tmax(0.0),
        xb_c(std::numeric_limits<double>::quiet_NaN()),
        q2_c(std::numeric_limits<double>::quiet_NaN()),
        t_c(std::numeric_limits<double>::quiet_NaN()) {}
};

struct FitResult {
    bool ok;
    double N;
    double N_stat;
    double phi_eval_deg;
    double xs_fit_phi_eval;
    double xs_fit_phi_eval_stat;
    double bh_phi_eval;
    double km15_phi_eval;
    double bh_over_km15_phi_eval;

    FitResult() :
        ok(false),
        N(0.0),
        N_stat(0.0),
        phi_eval_deg(0.001),
        xs_fit_phi_eval(0.0),
        xs_fit_phi_eval_stat(0.0),
        bh_phi_eval(0.0),
        km15_phi_eval(0.0),
        bh_over_km15_phi_eval(0.0) {}
};

static FitResult fit_norm_to_bh_and_eval_at_phi(const KinBinData &kb,
                                                double Ebeam,
                                                Helicity hel,
                                                double phi_eval_deg) {
    FitResult fr;
    fr.ok = false;
    fr.phi_eval_deg = phi_eval_deg;

    if (!finite_pos(kb.xb_c) || !finite_pos(kb.q2_c) || !finite_pos(kb.t_c)) {
        return fr;
    }

    // Weighted least squares for model: xs_i = N * BH_i
    //
    // Minimize sum_i (xs_i - N*BH_i)^2 / err_i^2
    //
    // N = sum_i (BH_i * xs_i / err_i^2) / sum_i (BH_i^2 / err_i^2)
    // Var(N) = 1 / sum_i (BH_i^2 / err_i^2)
    double num = 0.0;
    double den = 0.0;
    int n_used = 0;

    for (size_t i = 0; i < kb.points.size(); ++i) {
        const PhiPoint &p = kb.points[i];
        if (!std::isfinite(p.phi_deg)) continue;
        if (!(p.xs > 0.0) || !std::isfinite(p.xs)) continue;
        if (!(p.xs_stat > 0.0) || !std::isfinite(p.xs_stat)) continue;

        const double bh = vgg_bh_only(kb.xb_c, kb.q2_c, kb.t_c, p.phi_deg, Ebeam);
        if (!finite_pos(bh)) continue;

        const double w = 1.0 / (p.xs_stat * p.xs_stat);
        num += w * bh * p.xs;
        den += w * bh * bh;
        n_used += 1;
    }

    if (!(den > 0.0) || n_used < 2) {
        return fr;
    }

    fr.N = num / den;
    fr.N_stat = std::sqrt(1.0 / den);

    const double bh0 = vgg_bh_only(kb.xb_c, kb.q2_c, kb.t_c, phi_eval_deg, Ebeam);
    fr.bh_phi_eval = (std::isfinite(bh0) ? bh0 : 0.0);

    const double km0 = km15_xs(kb.xb_c, kb.q2_c, kb.t_c, phi_eval_deg, Ebeam, hel);
    fr.km15_phi_eval = (std::isfinite(km0) ? km0 : 0.0);

    if (finite_pos(fr.bh_phi_eval) && finite_pos(fr.km15_phi_eval)) {
        fr.bh_over_km15_phi_eval = fr.bh_phi_eval / fr.km15_phi_eval;
    } else {
        fr.bh_over_km15_phi_eval = 0.0;
    }

    if (finite_pos(fr.bh_phi_eval) && std::isfinite(fr.N) && fr.N > 0.0) {
        fr.xs_fit_phi_eval = fr.N * fr.bh_phi_eval;
        fr.xs_fit_phi_eval_stat = fr.N_stat * fr.bh_phi_eval;
    } else {
        fr.xs_fit_phi_eval = 0.0;
        fr.xs_fit_phi_eval_stat = 0.0;
    }

    fr.ok = std::isfinite(fr.N) && fr.N > 0.0 && std::isfinite(fr.N_stat) && fr.N_stat >= 0.0;
    return fr;
}

static void make_normalization_plots_legacy(const std::string &out_dir,
                                            const std::string &label,
                                            const std::string &helicity,
                                            const std::vector<PlotPoint> &pts) {
    ensure_output_dir_or_throw(out_dir);

    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);

    const std::string label_tag = sanitize_for_filename(label);
    const std::string hel_tag   = sanitize_for_filename(helicity);

    // Legacy behavior: keep x_min at 0.1 for log-x plots.
    const double x_min = 0.1;
    double x_max = x_min * 10.0;
    for (size_t i = 0; i < pts.size(); ++i) {
        if (std::isfinite(pts[i].d_edge_deg) && pts[i].d_edge_deg > x_max) {
            x_max = pts[i].d_edge_deg;
        }
    }
    x_max *= 1.10;

    // -------------------------------------------------------------------------
    // Legacy 1D: xs/BH vs d_edge, grouped by BH/KM15 bins
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
        } else {
            std::cerr << "[overall_norm] WARNING: no points available for legacy 1D plot.\n";
        }

        for (size_t si = 0; si < graphs.size(); ++si) {
            delete graphs[si];
        }
        delete c;
    }

    // -------------------------------------------------------------------------
    // Legacy 2D: BH/KM15 vs d_edge, colored by xs/BH
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
            std::cerr << "[overall_norm] WARNING: no points available for legacy 2D plot.\n";
        }

        delete g2;
    }
}

static std::vector<std::string> sorted_keys_by_numeric_minmax(const std::map<std::string, std::pair<double, double> > &key_to_edges) {
    std::vector<std::string> keys;
    keys.reserve(key_to_edges.size());
    for (std::map<std::string, std::pair<double, double> >::const_iterator it = key_to_edges.begin();
         it != key_to_edges.end(); ++it) {
        keys.push_back(it->first);
    }

    std::sort(keys.begin(), keys.end(),
              [&key_to_edges](const std::string &a, const std::string &b) {
                  const std::pair<double, double> ea = key_to_edges.find(a)->second;
                  const std::pair<double, double> eb = key_to_edges.find(b)->second;
                  if (ea.first < eb.first) return true;
                  if (ea.first > eb.first) return false;
                  return ea.second < eb.second;
              });

    return keys;
}

static void make_grid_canvases(const std::string &out_dir,
                               const std::string &label,
                               const std::string &helicity,
                               double Ebeam,
                               Helicity hel,
                               const std::map<std::string, KinBinData> &bins,
                               const std::map<std::string, FitResult> &fits) {
    ensure_output_dir_or_throw(out_dir);

    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);

    const std::string label_tag = sanitize_for_filename(label);
    const std::string hel_tag   = sanitize_for_filename(helicity);

    // Group kinematic bins by xB edges
    std::map<std::string, std::vector<std::string> > xb_to_bin_keys;
    std::map<std::string, std::pair<double, double> > xb_key_edges;

    for (std::map<std::string, KinBinData>::const_iterator it = bins.begin();
         it != bins.end(); ++it) {
        const std::string &kin_key = it->first;
        const KinBinData &kb = it->second;

        const std::string xb_key = xb_key_for_edges(kb.xbmin_s, kb.xbmax_s);
        xb_to_bin_keys[xb_key].push_back(kin_key);
        xb_key_edges[xb_key] = std::make_pair(kb.xbmin, kb.xbmax);
    }

    const std::vector<std::string> xb_keys_sorted = sorted_keys_by_numeric_minmax(xb_key_edges);

    for (size_t ix = 0; ix < xb_keys_sorted.size(); ++ix) {
        const std::string &xb_key = xb_keys_sorted[ix];
        const std::vector<std::string> &kin_keys = xb_to_bin_keys[xb_key];

        // For this xB slice, collect unique Q2 and |t| bins
        std::map<std::string, std::pair<double, double> > q2_key_edges;
        std::map<std::string, std::pair<double, double> > t_key_edges;

        // And build a lookup from (q2_key, t_key) -> kin_key
        std::map<std::string, std::string> q2t_to_kin;

        // Also keep representative edge strings for labeling
        std::map<std::string, std::pair<std::string, std::string> > q2_key_edge_strings;
        std::map<std::string, std::pair<std::string, std::string> > t_key_edge_strings;

        double xbmin_val = xb_key_edges[xb_key].first;
        double xbmax_val = xb_key_edges[xb_key].second;

        for (size_t k = 0; k < kin_keys.size(); ++k) {
            const std::string &kin_key = kin_keys[k];
            std::map<std::string, KinBinData>::const_iterator itb = bins.find(kin_key);
            if (itb == bins.end()) continue;

            const KinBinData &kb = itb->second;

            const std::string q2_key = q2_key_for_edges(kb.q2min_s, kb.q2max_s);
            const std::string t_key  = t_key_for_edges(kb.tmin_s, kb.tmax_s);

            q2_key_edges[q2_key] = std::make_pair(kb.q2min, kb.q2max);
            t_key_edges[t_key]   = std::make_pair(kb.tmin, kb.tmax);

            q2_key_edge_strings[q2_key] = std::make_pair(kb.q2min_s, kb.q2max_s);
            t_key_edge_strings[t_key]   = std::make_pair(kb.tmin_s,  kb.tmax_s);

            const std::string q2t_key = q2_key + "||" + t_key;
            q2t_to_kin[q2t_key] = kin_key;
        }

        const std::vector<std::string> q2_keys_sorted = sorted_keys_by_numeric_minmax(q2_key_edges);
        const std::vector<std::string> t_keys_sorted  = sorted_keys_by_numeric_minmax(t_key_edges);

        const int nrows = (int)q2_keys_sorted.size();
        const int ncols = (int)t_keys_sorted.size();

        if (nrows <= 0 || ncols <= 0) {
            std::cerr << "[overall_norm] WARNING: no (Q2,|t|) bins found for xB slice " << xb_key << "\n";
            continue;
        }

        // Canvas sizing: keep readable, similar scale per pad.
        const int pad_w = 360;
        const int pad_h = 300;
        const int W = pad_w * ncols + 200;
        const int H = pad_h * nrows + 180;

        const std::string cname = "c_grid_xb_" + sanitize_for_filename(xb_key) + "_" + label_tag + "_" + hel_tag;
        TCanvas *c = new TCanvas(cname.c_str(), cname.c_str(), W, H);
        c->SetLeftMargin(0.06);
        c->SetRightMargin(0.02);
        c->SetBottomMargin(0.06);
        c->SetTopMargin(0.08);

        c->Divide(ncols, nrows, 0.001, 0.001);

        // Title at top (NDC on canvas): we do it once using a transparent overlay pad.
        TPad *title_pad = new TPad("title_pad", "title_pad", 0.0, 0.92, 1.0, 1.0);
        title_pad->SetFillStyle(4000);
        title_pad->SetFrameFillStyle(4000);
        title_pad->SetLeftMargin(0.06);
        title_pad->SetRightMargin(0.02);
        title_pad->SetTopMargin(0.10);
        title_pad->SetBottomMargin(0.10);
        title_pad->Draw();
        title_pad->cd();

        TLatex ttl;
        ttl.SetNDC(kTRUE);
        ttl.SetTextSize(0.45);
        ttl.DrawLatex(0.02, 0.40,
                      Form("Normalization grids: %s, %s   xB in [%.3f, %.3f]   Ebeam=%.1f",
                           label.c_str(), helicity.c_str(), xbmin_val, xbmax_val, Ebeam));

        c->cd();

        // Build pads: row major with Q2 increasing downward is typical, but ROOT divide uses:
        // pad index = col + ncols*(row-1), with row starting at 1 from top.
        // We will map row=0 (lowest Q2) to bottom or top? Prefer top is lowest Q2 or vice versa.
        // For determinism, use Q2 sorted ascending with row=0 at top (ROOT's row 1 is top).
        for (int ir = 0; ir < nrows; ++ir) {
            const std::string &q2_key = q2_keys_sorted[(size_t)ir];
            const std::pair<std::string, std::string> q2_es = q2_key_edge_strings[q2_key];

            for (int ic = 0; ic < ncols; ++ic) {
                const std::string &t_key = t_keys_sorted[(size_t)ic];
                const std::pair<std::string, std::string> t_es = t_key_edge_strings[t_key];

                const int pad_idx = 1 + ic + ncols * ir;
                TPad *p = (TPad*)c->cd(pad_idx);
                if (!p) continue;

                p->SetLeftMargin(0.14);
                p->SetRightMargin(0.06);
                p->SetBottomMargin(0.16);
                p->SetTopMargin(0.18);
                p->SetTicks(1, 1);

                const std::string q2t_key = q2_key + "||" + t_key;
                std::map<std::string, std::string>::const_iterator itk = q2t_to_kin.find(q2t_key);

                // Empty pad if no bin exists
                if (itk == q2t_to_kin.end()) {
                    TLatex te;
                    te.SetNDC(kTRUE);
                    te.SetTextSize(0.070);
                    te.DrawLatex(0.15, 0.75, "No bin");
                    te.SetTextSize(0.055);
                    te.DrawLatex(0.15, 0.62, Form("Q2 in [%s,%s]", q2_es.first.c_str(), q2_es.second.c_str()));
                    te.DrawLatex(0.15, 0.50, Form("|t| in [%s,%s]", t_es.first.c_str(), t_es.second.c_str()));
                    continue;
                }

                const std::string &kin_key = itk->second;

                std::map<std::string, KinBinData>::const_iterator itb = bins.find(kin_key);
                std::map<std::string, FitResult>::const_iterator itf = fits.find(kin_key);

                if (itb == bins.end() || itf == fits.end() || !itf->second.ok) {
                    TLatex te;
                    te.SetNDC(kTRUE);
                    te.SetTextSize(0.070);
                    te.DrawLatex(0.15, 0.75, "Fit failed");
                    te.SetTextSize(0.055);
                    te.DrawLatex(0.15, 0.62, Form("Q2 in [%s,%s]", q2_es.first.c_str(), q2_es.second.c_str()));
                    te.DrawLatex(0.15, 0.50, Form("|t| in [%s,%s]", t_es.first.c_str(), t_es.second.c_str()));
                    continue;
                }

                const KinBinData &kb = itb->second;
                const FitResult &fr  = itf->second;

                // Build data graph (xs vs phi)
                TGraphErrors *gr = new TGraphErrors();
                gr->SetName(Form("gr_data_%d_%d", ir, ic));
                gr->SetMarkerStyle(20);
                gr->SetMarkerSize(0.8);
                gr->SetMarkerColor(kBlack);
                gr->SetLineColor(kBlack);
                gr->SetLineWidth(1);

                int np = 0;
                double y_max = 0.0;
                for (size_t ip = 0; ip < kb.points.size(); ++ip) {
                    const PhiPoint &pp = kb.points[ip];
                    if (!std::isfinite(pp.phi_deg)) continue;
                    if (!(pp.xs > 0.0) || !std::isfinite(pp.xs)) continue;
                    if (!(pp.xs_stat >= 0.0) || !std::isfinite(pp.xs_stat)) continue;

                    gr->SetPoint(np, pp.phi_deg, pp.xs);
                    gr->SetPointError(np, 0.0, pp.xs_stat);
                    if (pp.xs > y_max) y_max = pp.xs;
                    np += 1;
                }

                // Model curve: xs_model(phi) = N * BH(phi)
                // Use a modest sampling to keep it light but smooth.
                TGraph *gm = new TGraph();
                gm->SetName(Form("gr_model_%d_%d", ir, ic));
                gm->SetLineColor(kRed+1);
                gm->SetLineWidth(2);

                int nm = 0;
                for (int ph = 0; ph <= 360; ph += 2) {
                    const double phi_deg = (double)ph;
                    const double bh = vgg_bh_only(kb.xb_c, kb.q2_c, kb.t_c, phi_deg, Ebeam);
                    if (!finite_pos(bh)) continue;
                    const double y = fr.N * bh;
                    gm->SetPoint(nm, phi_deg, y);
                    if (y > y_max) y_max = y;
                    nm += 1;
                }

                // Axis frame via drawing data graph first
                gr->Draw("AP");
                gr->SetTitle("");
                if (gr->GetXaxis()) {
                    gr->GetXaxis()->SetTitle("#phi (deg)");
                    gr->GetXaxis()->SetLimits(0.0, 360.0);
                    gr->GetXaxis()->SetNdivisions(505);
                    gr->GetXaxis()->SetLabelSize(0.055);
                    gr->GetXaxis()->SetTitleSize(0.065);
                    gr->GetXaxis()->CenterTitle(true);
                }
                if (gr->GetYaxis()) {
                    gr->GetYaxis()->SetTitle("cross section");
                    gr->GetYaxis()->SetLabelSize(0.055);
                    gr->GetYaxis()->SetTitleSize(0.065);
                    gr->GetYaxis()->CenterTitle(true);
                    gr->GetYaxis()->SetRangeUser(0.0, (y_max > 0.0 ? 1.25 * y_max : 1.0));
                }

                gm->Draw("L same");
                gr->Draw("P same");

                // Mark phi evaluation point explicitly
                const double phi0 = fr.phi_eval_deg;
                TLine *lphi = new TLine(phi0, 0.0, phi0, (y_max > 0.0 ? 1.25 * y_max : 1.0));
                lphi->SetLineStyle(2);
                lphi->SetLineColor(kGray+2);
                lphi->SetLineWidth(2);
                lphi->Draw("same");

                // Annotate kinematics and fit output
                TLatex tx;
                tx.SetNDC(kTRUE);
                tx.SetTextSize(0.055);

                // Pad header
                tx.DrawLatex(0.14, 0.92, Form("Q2 in [%s,%s]   |t| in [%s,%s]", q2_es.first.c_str(), q2_es.second.c_str(),
                                              t_es.first.c_str(), t_es.second.c_str()));

                // Key output requested:
                // - Fitted normalization N +/- stat
                // - BH/KM15 evaluated at phi=0.001 deg
                tx.SetTextSize(0.060);
                tx.DrawLatex(0.14, 0.82, Form("Fit: xs(#phi) = N * BH(#phi)"));
                tx.DrawLatex(0.14, 0.74, Form("N = %.6f +/- %.6f", fr.N, fr.N_stat));
                tx.DrawLatex(0.14, 0.66, Form("BH/KM15 at #phi=%.3f deg: %.6f", fr.phi_eval_deg, fr.bh_over_km15_phi_eval));

                // Also show the fitted xs at phi=0.001 deg (often useful when eyeballing scale)
                tx.SetTextSize(0.055);
                tx.DrawLatex(0.14, 0.58, Form("xs_fit(#phi=%.3f deg) = %.3e +/- %.3e",
                                              fr.phi_eval_deg, fr.xs_fit_phi_eval, fr.xs_fit_phi_eval_stat));

                // Clean up per-pad heap objects (ROOT owns some drawables, but we delete explicitly)
                delete lphi;
                delete gm;
                delete gr;
            }
        }

        const std::string out_png = out_dir + "/norm_grid_xb_" + sanitize_for_filename(xb_key) + "_" + label_tag + "_" + hel_tag + ".png";
        c->SaveAs(out_png.c_str());

        delete title_pad;
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

        // Hard-coded evaluation point (deg), per your request.
        const double phi_eval_deg = 0.001;

        std::cout << "============================================================\n";
        std::cout << "[overall_norm] BH normalization study\n";
        std::cout << "[overall_norm] CSV      : " << csv_path << "\n";
        std::cout << "[overall_norm] label    : " << label << "\n";
        std::cout << "[overall_norm] helicity : " << helicity << "\n";
        std::cout << "[overall_norm] Ebeam    : " << Ebeam << "\n";
        std::cout << "[overall_norm] phi_eval : " << std::fixed << std::setprecision(3) << phi_eval_deg << " (deg)\n";
        std::cout << "------------------------------------------------------------\n";

        // Build full (xB,Q2,|t|) bins with all phi points for fitting.
        std::map<std::string, KinBinData> bins;

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

            KinBinData &kb = bins[key];
            if (kb.points.empty()) {
                kb.xbmin_s = xbmin_s;
                kb.xbmax_s = xbmax_s;
                kb.q2min_s = q2min_s;
                kb.q2max_s = q2max_s;
                kb.tmin_s  = tmin_s;
                kb.tmax_s  = tmax_s;

                kb.xbmin = std::atof(xbmin_s.c_str());
                kb.xbmax = std::atof(xbmax_s.c_str());
                kb.q2min = std::atof(q2min_s.c_str());
                kb.q2max = std::atof(q2max_s.c_str());
                kb.tmin  = std::atof(tmin_s.c_str());
                kb.tmax  = std::atof(tmax_s.c_str());

                kb.xb_c = xb_c;
                kb.q2_c = q2_c;
                kb.t_c  = t_c;
            }

            PhiPoint pp;
            pp.phi_deg = phi;
            pp.xs = xs.value;
            pp.xs_stat = xs.stat;
            pp.dist_to_edge = min_dist_to_0_or_360(phi);

            kb.points.push_back(pp);
        }

        if (bins.empty()) {
            std::cerr << "[overall_norm] WARNING: no bins found with positive cross section values.\n";
            return true;
        }

        // Derive legacy "best phi near 0/360" selection for the original 1D/2D plots and printout.
        std::map<std::string, BestPhiRow> best;
        for (std::map<std::string, KinBinData>::const_iterator it = bins.begin();
             it != bins.end(); ++it) {
            const std::string &key = it->first;
            const KinBinData &kb = it->second;

            BestPhiRow br;
            br.xb_c = kb.xb_c;
            br.q2_c = kb.q2_c;
            br.t_c  = kb.t_c;

            // Pick the phi point closest to 0 or 360 (true distance, no floor).
            for (size_t i = 0; i < kb.points.size(); ++i) {
                const PhiPoint &pp = kb.points[i];
                if (!std::isfinite(pp.phi_deg)) continue;
                if (!(pp.xs > 0.0) || !std::isfinite(pp.xs)) continue;

                if (pp.dist_to_edge < br.dist_to_edge) {
                    br.phi_deg = pp.phi_deg;
                    br.dist_to_edge = pp.dist_to_edge;
                    br.xs = pp.xs;
                    br.xs_stat = pp.xs_stat;
                    br.csv_row_index = 0; // not used in this version
                }
            }

            if (finite_pos(br.xb_c) && finite_pos(br.q2_c) && finite_pos(br.t_c) &&
                std::isfinite(br.phi_deg) && std::isfinite(br.dist_to_edge) &&
                br.xs > 0.0) {
                best[key] = br;
            }
        }

        std::cout << "[overall_norm] Unique kinematic bins found: " << bins.size()
                  << " (from " << n_rows << " CSV data rows)\n\n";

        // Legacy table header
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

        // Weighted mean accumulator for legacy summary (based on best-phi selection):
        // xs/BH for BH/KM15 in [0.95, 1.05] inclusive.
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

            // Legacy weighted mean summary (based on best-phi selection)
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

        // Legacy plots (unchanged behavior/ranges)
        const std::string out_dir_legacy = "output/normalization_study";
        make_normalization_plots_legacy(out_dir_legacy, label, helicity, plot_pts);

        // New: per-bin fit to xs(phi) = N * BH(phi) using all phi bins, and evaluation at phi=0.001 deg.
        std::map<std::string, FitResult> fits;
        fits.clear();

        for (std::map<std::string, KinBinData>::const_iterator it = bins.begin();
             it != bins.end(); ++it) {
            const std::string &kin_key = it->first;
            const KinBinData &kb = it->second;

            FitResult fr = fit_norm_to_bh_and_eval_at_phi(kb, Ebeam, hel, phi_eval_deg);
            fits[kin_key] = fr;
        }

        // Print a dedicated summary for the fitted normalization constants and BH/KM15 at phi=0.001 deg
        std::cout << "\n";
        std::cout << "------------------------------------------------------------\n";
        std::cout << "[overall_norm] Per-bin fits: xs(#phi) = N * BH(#phi)\n";
        std::cout << "[overall_norm] Evaluation at #phi = " << std::fixed << std::setprecision(3) << phi_eval_deg << " (deg)\n";
        std::cout
            << std::setw(8)  << "xB"
            << std::setw(10) << "Q2"
            << std::setw(12) << "-t"
            << std::setw(14) << "N"
            << std::setw(14) << "dN(stat)"
            << std::setw(16) << "BH/KM15(phi0)"
            << "\n";
        std::cout << std::string(8+10+12+14+14+16, '-') << "\n";

        for (std::map<std::string, KinBinData>::const_iterator it = bins.begin();
             it != bins.end(); ++it) {
            const std::string &kin_key = it->first;
            const KinBinData &kb = it->second;

            const FitResult &fr = fits[kin_key];
            if (!fr.ok) continue;

            const double tneg = -kb.t_c;

            std::cout
                << std::setw(8)  << std::fixed << std::setprecision(3) << kb.xb_c
                << std::setw(10) << std::fixed << std::setprecision(2) << kb.q2_c
                << std::setw(12) << std::fixed << std::setprecision(3) << tneg
                << std::setw(14) << std::fixed << std::setprecision(6) << fr.N
                << std::setw(14) << std::fixed << std::setprecision(6) << fr.N_stat
                << std::setw(16) << std::fixed << std::setprecision(6) << fr.bh_over_km15_phi_eval
                << "\n";
        }

        std::cout << "------------------------------------------------------------\n";

        // New grid canvases (xB canvases, rows=Q2, cols=|t|), each pad shows xs(phi), N*BH(phi),
        // and annotates N +/- stat and BH/KM15 at phi=0.001 deg.
        const std::string out_dir_grids = "output/normalization_study/grids";
        make_grid_canvases(out_dir_grids, label, helicity, Ebeam, hel, bins, fits);

        // Legacy weighted mean summary (unchanged)
        std::cout << "\n";
        std::cout << "------------------------------------------------------------\n";
        std::cout << "[overall_norm] Legacy weighted xs/BH for BH/KM15 in [0.95, 1.05]\n";
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