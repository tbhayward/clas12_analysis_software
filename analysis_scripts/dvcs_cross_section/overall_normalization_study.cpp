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
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <filesystem>

// ROOT
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TROOT.h>
#include <TPad.h>
#include <TH1F.h>
#include <TF1.h>

namespace {

// -----------------------------------------------------------------------------
// USER TOGGLE (edit here)
//   - true  : FAST PATH (one row per (xB,Q2,|t|) cell; pick phi closest to 0/360)
//   - false : USE ALL POINTS (every phi row in the CSV that has a valid xs)
// -----------------------------------------------------------------------------
static const bool kUseOnlyClosestToEdgePerKinCell = false;

// If kUseOnlyClosestToEdgePerKinCell == false, you can also optionally reject
// points with phi exactly on the edge (d_edge == 0). Keep as true by default.
static const bool kRequirePositiveDEdge = true;

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

static double min_dist_to_0_or_360(double phi_deg) {
    double d0   = std::fabs(phi_deg - 0.0);
    double d360 = std::fabs(phi_deg - 360.0);
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
// Grouping by BH/model
// -----------------------------------------------------------------------------

enum class BhGroupMetric {
    KM15 = 0,
    VGG  = 1
};

enum class BhKmGroup {
    G_95_105 = 0,
    G_90_95_OR_105_110,
    G_85_90_OR_110_115,
    G_OUTSIDE_15PCT,
    G_INVALID
};

static BhKmGroup categorize_bh_over_model(double r) {
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

static std::string ratio_name_for_metric(BhGroupMetric m) {
    if (m == BhGroupMetric::VGG) return "BH/VGG";
    return "BH/KM15";
}

static std::vector<GroupStyle> make_group_styles(const std::string &ratio_name) {
    std::vector<GroupStyle> s;
    s.push_back({BhKmGroup::G_95_105,             20, kBlack,     ratio_name + " in [0.95, 1.05]"});
    s.push_back({BhKmGroup::G_90_95_OR_105_110,   21, kBlue+1,    ratio_name + " in [0.90, 0.95) or (1.05, 1.10]"});
    s.push_back({BhKmGroup::G_85_90_OR_110_115,   22, kRed+1,     ratio_name + " in (0.85, 0.90) or (1.10, 1.15)"});
    s.push_back({BhKmGroup::G_OUTSIDE_15PCT,      24, kMagenta+1, ratio_name + " <= 0.85 or >= 1.15"});
    return s;
}

static int style_index_for_group(BhKmGroup g, const std::vector<GroupStyle> &styles) {
    for (size_t i = 0; i < styles.size(); ++i) {
        if (styles[i].group == g) return (int)i;
    }
    return -1;
}

// -----------------------------------------------------------------------------
// Data containers
// -----------------------------------------------------------------------------

struct PlotPoint {
    double d_edge_deg;
    double xs_over_bh;
    double xs_over_bh_err;
    double bh_over_km15;
    double bh_over_vgg;
};

struct DepPoint {
    double x;
    double y;
    double ey;
    double bh_over_km15;
    double bh_over_vgg;
};

static double ratio_for_metric(const PlotPoint &p, BhGroupMetric m) {
    return (m == BhGroupMetric::KM15 ? p.bh_over_km15 : p.bh_over_vgg);
}

static double ratio_for_metric(const DepPoint &p, BhGroupMetric m) {
    return (m == BhGroupMetric::KM15 ? p.bh_over_km15 : p.bh_over_vgg);
}

// -----------------------------------------------------------------------------
// d_edge plot: 1x2
// -----------------------------------------------------------------------------

static void draw_dedge_panel(TPad *pad,
                            const std::vector<PlotPoint> &pts,
                            const std::string &label,
                            const std::string &helicity,
                            double x_min,
                            double x_max,
                            BhGroupMetric metric) {
    pad->SetLogx(1);
    pad->SetLeftMargin(0.14);
    pad->SetRightMargin(0.05);
    pad->SetBottomMargin(0.12);
    pad->SetTopMargin(0.12);
    pad->SetGrid(1, 1);

    const std::string ratio_name = ratio_name_for_metric(metric);
    const std::vector<GroupStyle> styles = make_group_styles(ratio_name);

    TH1F *frame = (TH1F*)pad->DrawFrame(x_min, 0.0, x_max, 2.0);
    frame->SetTitle("");
    frame->GetXaxis()->SetTitle("d_edge (deg)");
    frame->GetYaxis()->SetTitle("xs/BH");
    frame->GetXaxis()->SetNdivisions(505);
    frame->GetXaxis()->CenterTitle();
    frame->GetYaxis()->CenterTitle();

    TLatex tl;
    tl.SetNDC(kTRUE);
    tl.SetTextSize(0.030);
    tl.DrawLatex(0.14, 0.93, Form("%s, %s   (%s)", label.c_str(), helicity.c_str(), ratio_name.c_str()));

    TLine *l = new TLine(x_min, 1.0, x_max, 1.0);
    l->SetLineStyle(2);
    l->SetLineColor(kGray+2);
    l->SetLineWidth(2);
    l->Draw("same");
    l->SetBit(kCanDelete);

    std::vector<TGraphErrors*> graphs(styles.size(), (TGraphErrors*)nullptr);
    std::vector<int> n_in(styles.size(), 0);

    for (size_t si = 0; si < styles.size(); ++si) {
        graphs[si] = new TGraphErrors();
        graphs[si]->SetName(Form("gr_dedge_%s_%zu", sanitize_for_filename(ratio_name).c_str(), si));
        graphs[si]->SetMarkerStyle(styles[si].marker_style);
        graphs[si]->SetMarkerColor(styles[si].marker_color);
        graphs[si]->SetLineColor(styles[si].marker_color);
        graphs[si]->SetLineWidth(1);
        graphs[si]->SetMarkerSize(0.85);
    } //endfor

    for (size_t i = 0; i < pts.size(); ++i) {
        const double r = ratio_for_metric(pts[i], metric);
        const BhKmGroup g = categorize_bh_over_model(r);
        const int si = style_index_for_group(g, styles);
        if (si < 0) continue;

        const int n = n_in[(size_t)si];
        graphs[(size_t)si]->SetPoint(n, pts[i].d_edge_deg, pts[i].xs_over_bh);
        graphs[(size_t)si]->SetPointError(n, 0.0, pts[i].xs_over_bh_err);
        n_in[(size_t)si] += 1;
    } //endfor

    bool any = false;
    for (size_t si = 0; si < graphs.size(); ++si) {
        if (graphs[si] && graphs[si]->GetN() > 0) {
            graphs[si]->Draw("PE1 SAME");
            graphs[si]->SetBit(kCanDelete);
            any = true;
        } else {
            delete graphs[si];
            graphs[si] = nullptr;
        }
    } //endfor

    if (!any) {
        TLatex t0;
        t0.SetNDC(kTRUE);
        t0.SetTextSize(0.050);
        t0.DrawLatex(0.18, 0.70, "No points");
        pad->Update();
        return;
    } //endif

    TLegend *leg = new TLegend(0.16, 0.62, 0.74, 0.88);
    leg->SetBorderSize(1);
    leg->SetFillStyle(1001);
    leg->SetFillColor(kWhite);
    leg->SetTextSize(0.024);

    for (size_t si = 0; si < graphs.size(); ++si) {
        if (graphs[si] && graphs[si]->GetN() > 0) {
            leg->AddEntry(graphs[si], styles[si].label.c_str(), "p");
        }
    } //endfor
    leg->Draw();
    leg->SetBit(kCanDelete);

    pad->Update();
}

static void make_normalization_plot_dedge_1x2(const std::string &out_dir,
                                              const std::string &label,
                                              const std::string &helicity,
                                              const std::vector<PlotPoint> &pts) {
    ensure_output_dir_or_throw(out_dir);

    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    const std::string label_tag = sanitize_for_filename(label);
    const std::string hel_tag   = sanitize_for_filename(helicity);

    const double x_min = 1.0;
    double x_max = x_min * 2.0;

    for (size_t i = 0; i < pts.size(); ++i) {
        if (std::isfinite(pts[i].d_edge_deg) && pts[i].d_edge_deg > x_max) {
            x_max = pts[i].d_edge_deg;
        }
    }
    x_max *= 1.10;

    TCanvas *c = new TCanvas("c_norm_dedge_1x2", "c_norm_dedge_1x2", 1350, 620);
    c->Divide(2, 1, 0.001, 0.001);

    {
        TPad *p = (TPad*)c->cd(1);
        draw_dedge_panel(p, pts, label, helicity, x_min, x_max, BhGroupMetric::KM15);
    }

    {
        TPad *p = (TPad*)c->cd(2);
        draw_dedge_panel(p, pts, label, helicity, x_min, x_max, BhGroupMetric::VGG);
    }

    const std::string out_png = out_dir + "/norm_1d_xs_over_bh_vs_dedge_" + label_tag + "_" + hel_tag + ".png";
    c->SaveAs(out_png.c_str());

    delete c;
}

// -----------------------------------------------------------------------------
// Dependence plots: 2x2
// -----------------------------------------------------------------------------

static void draw_all_groups_dep_pad(TPad *p,
                                    const std::vector<DepPoint> &pts,
                                    const std::string &label,
                                    const std::string &helicity,
                                    const std::string &x_title,
                                    double x_min,
                                    double x_max,
                                    BhGroupMetric metric) {
    p->SetLeftMargin(0.13);
    p->SetRightMargin(0.05);
    p->SetBottomMargin(0.14);
    p->SetTopMargin(0.12);
    p->SetGrid(1, 1);

    const std::string ratio_name = ratio_name_for_metric(metric);
    const std::vector<GroupStyle> styles = make_group_styles(ratio_name);

    TH1F *frame = (TH1F*)p->DrawFrame(x_min, 0.0, x_max, 2.0);
    frame->SetTitle("");
    frame->GetXaxis()->SetTitle(x_title.c_str());
    frame->GetYaxis()->SetTitle("xs/BH");
    frame->GetXaxis()->SetNdivisions(505);
    frame->GetXaxis()->CenterTitle();
    frame->GetYaxis()->CenterTitle();

    TLatex tl;
    tl.SetNDC(kTRUE);
    tl.SetTextSize(0.030);
    tl.DrawLatex(0.12, 0.93, Form("%s, %s   (%s)", label.c_str(), helicity.c_str(), ratio_name.c_str()));

    std::vector<TGraphErrors*> graphs(styles.size(), (TGraphErrors*)nullptr);
    std::vector<int> n_in(styles.size(), 0);

    for (size_t si = 0; si < styles.size(); ++si) {
        graphs[si] = new TGraphErrors();
        graphs[si]->SetName(Form("gr_dep_%s_%zu", sanitize_for_filename(ratio_name).c_str(), si));
        graphs[si]->SetMarkerStyle(styles[si].marker_style);
        graphs[si]->SetMarkerColor(styles[si].marker_color);
        graphs[si]->SetLineColor(styles[si].marker_color);
        graphs[si]->SetLineWidth(1);
        graphs[si]->SetMarkerSize(0.85);
    } //endfor

    for (size_t i = 0; i < pts.size(); ++i) {
        const double r = ratio_for_metric(pts[i], metric);
        const BhKmGroup g = categorize_bh_over_model(r);
        const int si = style_index_for_group(g, styles);
        if (si < 0) continue;

        const int n = n_in[(size_t)si];
        graphs[(size_t)si]->SetPoint(n, pts[i].x, pts[i].y);
        graphs[(size_t)si]->SetPointError(n, 0.0, pts[i].ey);
        n_in[(size_t)si] += 1;
    } //endfor

    bool any = false;
    for (size_t si = 0; si < graphs.size(); ++si) {
        if (graphs[si] && graphs[si]->GetN() > 0) {
            graphs[si]->Draw("PE1 SAME");
            graphs[si]->SetBit(kCanDelete);
            any = true;
        } else {
            delete graphs[si];
            graphs[si] = nullptr;
        }
    } //endfor

    if (!any) {
        TLatex t0;
        t0.SetNDC(kTRUE);
        t0.SetTextSize(0.050);
        t0.DrawLatex(0.18, 0.70, "No points");
        p->Update();
        return;
    } //endif

    TLegend *leg = new TLegend(0.16, 0.62, 0.78, 0.88);
    leg->SetBorderSize(1);
    leg->SetFillStyle(1001);
    leg->SetFillColor(kWhite);
    leg->SetTextSize(0.024);

    for (size_t si = 0; si < graphs.size(); ++si) {
        if (graphs[si] && graphs[si]->GetN() > 0) {
            leg->AddEntry(graphs[si], styles[si].label.c_str(), "p");
        }
    } //endfor
    leg->Draw();
    leg->SetBit(kCanDelete);

    p->Update();
}

static void draw_clean_fit_dep_pad(TPad *p,
                                   const std::vector<DepPoint> &pts,
                                   const std::string &label,
                                   const std::string &helicity,
                                   const std::string &x_title,
                                   double x_min,
                                   double x_max,
                                   BhGroupMetric metric) {
    p->SetLeftMargin(0.13);
    p->SetRightMargin(0.05);
    p->SetBottomMargin(0.14);
    p->SetTopMargin(0.12);
    p->SetGrid(1, 1);

    const std::string ratio_name = ratio_name_for_metric(metric);

    TH1F *frame = (TH1F*)p->DrawFrame(x_min, 0.0, x_max, 2.0);
    frame->SetTitle("");
    frame->GetXaxis()->SetTitle(x_title.c_str());
    frame->GetYaxis()->SetTitle("xs/BH");
    frame->GetXaxis()->SetNdivisions(505);
    frame->GetXaxis()->CenterTitle();
    frame->GetYaxis()->CenterTitle();

    TLatex tl;
    tl.SetNDC(kTRUE);
    tl.SetTextSize(0.030);
    tl.DrawLatex(0.12, 0.93, Form("%s, %s   (%s)", label.c_str(), helicity.c_str(), ratio_name.c_str()));

    TGraphErrors *gr_clean = new TGraphErrors();
    gr_clean->SetName(Form("gr_clean_%s", sanitize_for_filename(ratio_name).c_str()));
    gr_clean->SetMarkerStyle(20);
    gr_clean->SetMarkerColor(kBlack);
    gr_clean->SetLineColor(kBlack);
    gr_clean->SetLineWidth(1);
    gr_clean->SetMarkerSize(0.85);

    int n = 0;
    for (size_t i = 0; i < pts.size(); ++i) {
        const double r = ratio_for_metric(pts[i], metric);
        const BhKmGroup g = categorize_bh_over_model(r);
        if (g != BhKmGroup::G_95_105) continue;

        gr_clean->SetPoint(n, pts[i].x, pts[i].y);
        gr_clean->SetPointError(n, 0.0, pts[i].ey);
        n += 1;
    } //endfor

    if (gr_clean->GetN() <= 0) {
        delete gr_clean;
        TLatex t0;
        t0.SetNDC(kTRUE);
        t0.SetTextSize(0.050);
        t0.DrawLatex(0.18, 0.70, "No points in [0.95, 1.05]");
        p->Update();
        return;
    } //endif

    gr_clean->Draw("PE1 SAME");
    gr_clean->SetBit(kCanDelete);

    if (gr_clean->GetN() >= 2) {
        TF1 *f_lin = new TF1(Form("f_lin_%s", sanitize_for_filename(ratio_name).c_str()), "pol1", x_min, x_max);

        gr_clean->Fit(f_lin, "Q0");

        f_lin->SetLineColor(kRed+1);
        f_lin->SetLineStyle(2);
        f_lin->SetLineWidth(2);

        f_lin->Draw("SAME");
        f_lin->SetBit(kCanDelete);

        const double p0 = f_lin->GetParameter(0);
        const double p1 = f_lin->GetParameter(1);
        const double e0 = f_lin->GetParError(0);
        const double e1 = f_lin->GetParError(1);

        TLegend *legfit = new TLegend(0.58, 0.72, 0.93, 0.88);
        legfit->SetBorderSize(1);
        legfit->SetFillStyle(1001);
        legfit->SetFillColor(kWhite);
        legfit->SetTextSize(0.028);

        legfit->AddEntry((TObject*)0, "Linear fit (pol1)", "");
        legfit->AddEntry((TObject*)0, Form("p0 = %.4f #pm %.4f", p0, e0), "");
        legfit->AddEntry((TObject*)0, Form("p1 = %.4f #pm %.4f", p1, e1), "");

        legfit->Draw();
        legfit->SetBit(kCanDelete);
    } else {
        TLatex t1;
        t1.SetNDC(kTRUE);
        t1.SetTextSize(0.050);
        t1.DrawLatex(0.18, 0.70, "Too few points to fit");
    } //endif

    p->Update();
}

static void make_dependence_plot_2x2(const std::string &out_dir,
                                     const std::string &label,
                                     const std::string &helicity,
                                     const std::string &x_title,
                                     const std::string &file_tag,
                                     double x_min,
                                     double x_max,
                                     const std::vector<DepPoint> &pts) {
    ensure_output_dir_or_throw(out_dir);

    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    const std::string label_tag = sanitize_for_filename(label);
    const std::string hel_tag   = sanitize_for_filename(helicity);

    TCanvas *c = new TCanvas(Form("c_dep_%s_2x2", file_tag.c_str()),
                             Form("c_dep_%s_2x2", file_tag.c_str()),
                             1200, 980);
    c->Divide(2, 2, 0.001, 0.001);

    {
        TPad *p = (TPad*)c->cd(1);
        draw_all_groups_dep_pad(p, pts, label, helicity, x_title, x_min, x_max, BhGroupMetric::KM15);
    }

    {
        TPad *p = (TPad*)c->cd(2);
        draw_clean_fit_dep_pad(p, pts, label, helicity, x_title, x_min, x_max, BhGroupMetric::KM15);
    }

    {
        TPad *p = (TPad*)c->cd(3);
        draw_all_groups_dep_pad(p, pts, label, helicity, x_title, x_min, x_max, BhGroupMetric::VGG);
    }

    {
        TPad *p = (TPad*)c->cd(4);
        draw_clean_fit_dep_pad(p, pts, label, helicity, x_title, x_min, x_max, BhGroupMetric::VGG);
    }

    const std::string out_png = out_dir + "/norm_xs_over_bh_vs_" + file_tag + "_" + label_tag + "_" + hel_tag + ".png";
    c->SaveAs(out_png.c_str());

    delete c;
}

// -----------------------------------------------------------------------------
// Row container
//   - When using ALL points: we keep per-row values and compute BH/model per row.
//   - When using BEST points: we keep one per (xB,Q2,|t|) cell (closest phi edge).
// -----------------------------------------------------------------------------

struct RowData {
    double xb_c;
    double q2_c;
    double t_c;           // positive |t|
    double phi_deg;
    double dist_to_edge;  // min(|phi-0|, |phi-360|)
    double xs;
    double xs_stat;
    double e_theta;
    double p_theta;
    double g_theta;
    size_t csv_row_index;

    RowData() :
        xb_c(std::numeric_limits<double>::quiet_NaN()),
        q2_c(std::numeric_limits<double>::quiet_NaN()),
        t_c(std::numeric_limits<double>::quiet_NaN()),
        phi_deg(std::numeric_limits<double>::quiet_NaN()),
        dist_to_edge(std::numeric_limits<double>::infinity()),
        xs(0.0),
        xs_stat(0.0),
        e_theta(std::numeric_limits<double>::quiet_NaN()),
        p_theta(std::numeric_limits<double>::quiet_NaN()),
        g_theta(std::numeric_limits<double>::quiet_NaN()),
        csv_row_index(0) {}
};

// -----------------------------------------------------------------------------
// Main entry
// -----------------------------------------------------------------------------

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

        // FIX: theta columns do NOT have "_avg" in the CSV schema
        const std::string col_e_theta = "e_theta, " + label;
        const std::string col_p_theta = "p_theta, " + label;
        const std::string col_g_theta = "g_theta, " + label;

        const int c_e_theta = require_col(header, col_e_theta);
        const int c_p_theta = require_col(header, col_p_theta);
        const int c_g_theta = require_col(header, col_g_theta);

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
        std::cout << "[overall_norm] mode     : "
                  << (kUseOnlyClosestToEdgePerKinCell ? "closest-to-edge (per kinematic cell)" : "all points")
                  << "\n";
        std::cout << "------------------------------------------------------------\n";

        // ---------------------------------------------------------------------
        // Build the selected row set based on the toggle
        // ---------------------------------------------------------------------

        std::vector<RowData> selected_rows;
        selected_rows.reserve(lines.size() > 1 ? lines.size() - 1 : 0);

        if (kUseOnlyClosestToEdgePerKinCell) {
            // Map: key = (xBmin,xBmax,Q2min,Q2max,tmin,tmax), store best (smallest dist)
            std::map<std::string, RowData> best;

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

                const double xb_c = std::atof(trim(unquote(fields[c_xbavg])).c_str());
                const double q2_c = std::atof(trim(unquote(fields[c_q2avg])).c_str());
                const double t_c  = std::atof(trim(unquote(fields[c_tavg])).c_str());
                const double phi  = std::atof(trim(unquote(fields[c_phiavg])).c_str());

                const double e_th = std::atof(trim(unquote(fields[c_e_theta])).c_str());
                const double p_th = std::atof(trim(unquote(fields[c_p_theta])).c_str());
                const double g_th = std::atof(trim(unquote(fields[c_g_theta])).c_str());

                if (!finite_pos(xb_c) || !finite_pos(q2_c) || !finite_pos(t_c) || !std::isfinite(phi)) {
                    continue;
                } //endif

                TripleCell xs = parse_tuple3(fields[c_xs]);
                if (!(xs.value > 0.0) || !std::isfinite(xs.value)) {
                    continue;
                } //endif

                const double dist = min_dist_to_0_or_360(phi);
                if (kRequirePositiveDEdge && !(dist > 0.0)) {
                    continue;
                } //endif

                std::map<std::string, RowData>::iterator it = best.find(key);
                if (it == best.end() || dist < it->second.dist_to_edge) {
                    RowData rd;
                    rd.xb_c = xb_c;
                    rd.q2_c = q2_c;
                    rd.t_c  = t_c;
                    rd.phi_deg = phi;
                    rd.dist_to_edge = dist;
                    rd.xs = xs.value;
                    rd.xs_stat = xs.stat;
                    rd.e_theta = e_th;
                    rd.p_theta = p_th;
                    rd.g_theta = g_th;
                    rd.csv_row_index = r;
                    best[key] = rd;
                } //endif
            } //endfor rows

            for (std::map<std::string, RowData>::const_iterator it = best.begin(); it != best.end(); ++it) {
                selected_rows.push_back(it->second);
            } //endfor

            if (selected_rows.empty()) {
                std::cerr << "[overall_norm] WARNING: no selected bins found in closest-to-edge mode.\n";
                return true;
            } //endif

            std::cout << "[overall_norm] Selected rows: " << selected_rows.size()
                      << " (one per kinematic cell, closest to phi edge)\n\n";

        } else {
            // Use every valid row
            for (size_t r = 1; r < lines.size(); ++r) {
                if (lines[r].empty()) continue;

                std::vector<std::string> fields = split_csv_line(lines[r]);
                if (fields.size() != header.size()) continue;

                const double xb_c = std::atof(trim(unquote(fields[c_xbavg])).c_str());
                const double q2_c = std::atof(trim(unquote(fields[c_q2avg])).c_str());
                const double t_c  = std::atof(trim(unquote(fields[c_tavg])).c_str());
                const double phi  = std::atof(trim(unquote(fields[c_phiavg])).c_str());

                const double e_th = std::atof(trim(unquote(fields[c_e_theta])).c_str());
                const double p_th = std::atof(trim(unquote(fields[c_p_theta])).c_str());
                const double g_th = std::atof(trim(unquote(fields[c_g_theta])).c_str());

                if (!finite_pos(xb_c) || !finite_pos(q2_c) || !finite_pos(t_c) || !std::isfinite(phi)) {
                    continue;
                } //endif

                TripleCell xs = parse_tuple3(fields[c_xs]);
                if (!(xs.value > 0.0) || !std::isfinite(xs.value)) {
                    continue;
                } //endif

                const double dist = min_dist_to_0_or_360(phi);
                if (kRequirePositiveDEdge && !(dist > 0.0)) {
                    continue;
                } //endif

                RowData rd;
                rd.xb_c = xb_c;
                rd.q2_c = q2_c;
                rd.t_c  = t_c;
                rd.phi_deg = phi;
                rd.dist_to_edge = dist;
                rd.xs = xs.value;
                rd.xs_stat = xs.stat;
                rd.e_theta = e_th;
                rd.p_theta = p_th;
                rd.g_theta = g_th;
                rd.csv_row_index = r;

                selected_rows.push_back(rd);
            } //endfor rows

            if (selected_rows.empty()) {
                std::cerr << "[overall_norm] WARNING: no selected bins found in all-points mode.\n";
                return true;
            } //endif

            std::cout << "[overall_norm] Selected rows: " << selected_rows.size()
                      << " (all points)\n\n";
        } //endif mode switch

        // ---------------------------------------------------------------------
        // Print header for the table (same columns as before)
        // ---------------------------------------------------------------------
        std::cout
            << std::setw(8)  << "xB"
            << std::setw(10) << "Q2"
            << std::setw(12) << "t_abs"
            << std::setw(10) << "d_edge"
            << std::setw(14) << "xs"
            << std::setw(14) << "xs/BH"
            << std::setw(14) << "BH/VGG"
            << std::setw(14) << "BH/KM15"
            << "\n";
        std::cout << std::string(8+10+12+10+14*4, '-') << "\n";

        std::vector<PlotPoint> plot_pts_dedge;
        plot_pts_dedge.reserve(selected_rows.size());

        std::vector<DepPoint> dep_xb;
        std::vector<DepPoint> dep_q2;
        std::vector<DepPoint> dep_t;
        std::vector<DepPoint> dep_e_theta;
        std::vector<DepPoint> dep_p_theta;
        std::vector<DepPoint> dep_g_theta;

        dep_xb.reserve(selected_rows.size());
        dep_q2.reserve(selected_rows.size());
        dep_t.reserve(selected_rows.size());
        dep_e_theta.reserve(selected_rows.size());
        dep_p_theta.reserve(selected_rows.size());
        dep_g_theta.reserve(selected_rows.size());

        double sumw = 0.0;
        double sumwx = 0.0;
        int n_weighted_used = 0;
        int n_in_95_105_total = 0;

        for (size_t irow = 0; irow < selected_rows.size(); ++irow) {
            const RowData &br = selected_rows[irow];

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
                } //endif
            } //endif
            if (finite_pos(vgg)) {
                bh_over_vgg = bh / vgg;
            } //endif
            if (finite_pos(km)) {
                bh_over_km = bh / km;
            } //endif

            std::cout
                << std::setw(8)  << std::fixed << std::setprecision(3) << xb
                << std::setw(10) << std::fixed << std::setprecision(2) << q2
                << std::setw(12) << std::fixed << std::setprecision(3) << tpos
                << std::setw(10) << std::fixed << std::setprecision(1) << br.dist_to_edge
                << std::setw(14) << std::scientific << std::setprecision(3) << br.xs
                << std::setw(14) << std::fixed << std::setprecision(3) << xs_over_bh
                << std::setw(14) << std::fixed << std::setprecision(3) << bh_over_vgg
                << std::setw(14) << std::fixed << std::setprecision(3) << bh_over_km
                << "\n";

            // d_edge plot points: require dist_to_edge > 0.0 unless user disabled
            if ((!kRequirePositiveDEdge || (std::isfinite(br.dist_to_edge) && br.dist_to_edge > 0.0)) &&
                std::isfinite(xs_over_bh) && xs_over_bh >= 0.0 &&
                std::isfinite(bh_over_km) && bh_over_km > 0.0 &&
                std::isfinite(bh_over_vgg) && bh_over_vgg > 0.0) {

                PlotPoint p;
                p.d_edge_deg      = br.dist_to_edge;
                p.xs_over_bh      = xs_over_bh;
                p.xs_over_bh_err  = xs_over_bh_err;
                p.bh_over_km15    = bh_over_km;
                p.bh_over_vgg     = bh_over_vgg;
                plot_pts_dedge.push_back(p);
            } //endif

            // Dependence plots
            if (std::isfinite(xs_over_bh) && xs_over_bh >= 0.0 &&
                std::isfinite(bh_over_km) && bh_over_km > 0.0 &&
                std::isfinite(bh_over_vgg) && bh_over_vgg > 0.0) {

                DepPoint px;
                px.x = xb;
                px.y = xs_over_bh;
                px.ey = (std::isfinite(xs_over_bh_err) ? xs_over_bh_err : 0.0);
                px.bh_over_km15 = bh_over_km;
                px.bh_over_vgg  = bh_over_vgg;
                dep_xb.push_back(px);

                DepPoint pq;
                pq.x = q2;
                pq.y = xs_over_bh;
                pq.ey = (std::isfinite(xs_over_bh_err) ? xs_over_bh_err : 0.0);
                pq.bh_over_km15 = bh_over_km;
                pq.bh_over_vgg  = bh_over_vgg;
                dep_q2.push_back(pq);

                DepPoint pt;
                pt.x = tpos;
                pt.y = xs_over_bh;
                pt.ey = (std::isfinite(xs_over_bh_err) ? xs_over_bh_err : 0.0);
                pt.bh_over_km15 = bh_over_km;
                pt.bh_over_vgg  = bh_over_vgg;
                dep_t.push_back(pt);

                if (std::isfinite(br.e_theta)) {
                    DepPoint pe;
                    pe.x = br.e_theta;
                    pe.y = xs_over_bh;
                    pe.ey = (std::isfinite(xs_over_bh_err) ? xs_over_bh_err : 0.0);
                    pe.bh_over_km15 = bh_over_km;
                    pe.bh_over_vgg  = bh_over_vgg;
                    dep_e_theta.push_back(pe);
                } //endif

                if (std::isfinite(br.p_theta)) {
                    DepPoint pp;
                    pp.x = br.p_theta;
                    pp.y = xs_over_bh;
                    pp.ey = (std::isfinite(xs_over_bh_err) ? xs_over_bh_err : 0.0);
                    pp.bh_over_km15 = bh_over_km;
                    pp.bh_over_vgg  = bh_over_vgg;
                    dep_p_theta.push_back(pp);
                } //endif

                if (std::isfinite(br.g_theta)) {
                    DepPoint pg;
                    pg.x = br.g_theta;
                    pg.y = xs_over_bh;
                    pg.ey = (std::isfinite(xs_over_bh_err) ? xs_over_bh_err : 0.0);
                    pg.bh_over_km15 = bh_over_km;
                    pg.bh_over_vgg  = bh_over_vgg;
                    dep_g_theta.push_back(pg);
                } //endif
            } //endif

            if (std::isfinite(bh_over_km) && bh_over_km >= 0.95 && bh_over_km <= 1.05) {
                n_in_95_105_total += 1;

                if (std::isfinite(xs_over_bh_err) && xs_over_bh_err > 0.0 &&
                    std::isfinite(xs_over_bh)) {

                    const double w = 1.0 / (xs_over_bh_err * xs_over_bh_err);
                    sumw  += w;
                    sumwx += w * xs_over_bh;
                    n_weighted_used += 1;
                } //endif
            } //endif
        } //endfor selected rows

        // ---------------------------------------------------------------------
        // OUTPUT DIRECTORY POLICY:
        //   output/normalization_study/<RunPeriodLabelSanitized>
        // ---------------------------------------------------------------------
        const std::string out_root = "output/normalization_study";
        ensure_output_dir_or_throw(out_root);

        const std::string out_dir = out_root + "/" + sanitize_for_filename(label);
        ensure_output_dir_or_throw(out_dir);

        make_normalization_plot_dedge_1x2(out_dir, label, helicity, plot_pts_dedge);

        make_dependence_plot_2x2(out_dir, label, helicity,
                                 "xB",
                                 "xb",
                                 0.0, 0.6,
                                 dep_xb);

        make_dependence_plot_2x2(out_dir, label, helicity,
                                 "Q^{2}",
                                 "q2",
                                 1.0, 6.0,
                                 dep_q2);

        // Label changed from t_abs to -t (values are already positive)
        make_dependence_plot_2x2(out_dir, label, helicity,
                                 "-t",
                                 "t",
                                 0.0, 1.0,
                                 dep_t);

        make_dependence_plot_2x2(out_dir, label, helicity,
                                 "#theta_{e} (deg)",
                                 "e_theta",
                                 0.0, 40.0,
                                 dep_e_theta);

        make_dependence_plot_2x2(out_dir, label, helicity,
                                 "#theta_{p} (deg)",
                                 "p_theta",
                                 0.0, 80.0,
                                 dep_p_theta);

        make_dependence_plot_2x2(out_dir, label, helicity,
                                 "#theta_{#gamma} (deg)",
                                 "g_theta",
                                 0.0, 40.0,
                                 dep_g_theta);

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
        } //endif
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