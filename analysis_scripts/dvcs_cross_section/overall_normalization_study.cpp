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

static bool has_helicity_resolved_cross_sections(const std::string &label) {
    if (label == "Sp18 Inb" || label == "Sp18 Out" ||
        label == "Sp18" || label == "10.6 GeV") {
        return false;
    }

    return true;
}

static bool normalization_helicity_allowed(const std::string &label,
                                           const std::string &helicity) {
    if (helicity == "unpol") {
        return true;
    }

    return has_helicity_resolved_cross_sections(label);
}


// -----------------------------------------------------------------------------
// USER TOGGLES (edit here)
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// If true, skip ALL calculations and ALL plotting, and just write "1.00"
// into every row of "norm, <label>" if the column exists, then return.
// -----------------------------------------------------------------------------
static bool kSkipAllWorkWriteUnityNorm = true;

// -----------------------------------------------------------------------------
// Edge-window normalization sample.
//
// Default behavior:
//   - use all valid points with d_edge <= 10 deg
//   - do not use BH/VGG or BH/KM15 to decide which points enter the sample
//
// Here
//
//   d_edge = min(distance to phi = 0 deg, distance to phi = 360 deg)
//
// If kRequirePositiveDEdge is true, points exactly at d_edge = 0 are rejected.
// -----------------------------------------------------------------------------
static bool kUseAllPointsWithinEdgeWindow = true;
static bool kRequirePositiveDEdge = true;
static double kMaxDEdgeForNormalizationDeg = 10.0;

// -----------------------------------------------------------------------------
// Optional alternate mode.
//
// If kUseAllPointsWithinEdgeWindow is false, the code uses at most one point per
// (xB,Q2,|t|) cell, but still only among points that pass the edge-window cut.
// In other words, even the closest-to-edge mode will never use d_edge > 10 deg.
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// Toggle for which variable is used to parameterize the normalization written to
// "norm, <label>".
//
//   - PTheta : norm = p0 + p1 * p_theta
//   - XB     : norm = p0 + p1 * xBavg
//
// This only affects the "norm write" section at the end. All dependence plots
// are still produced.
// -----------------------------------------------------------------------------
enum class NormXAxis {
    PTheta = 0,
    XB     = 1
};

// Default: use xB to parameterize the normalization.
static NormXAxis kNormXAxis = NormXAxis::XB;
static std::string kOutputRoot = "output/normalization_study";

// -----------------------------------------------------------------------------
// Model-ratio diagnostic grouping.
//
// These ratios are printed and used only to style diagnostic plots. They are
// NOT used to select points for the normalization fit.
// -----------------------------------------------------------------------------
enum class BhGroupMetric {
    KM15 = 0,
    VGG  = 1
};

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
    } //endfor

    out.push_back(field);
    return out;
}

static std::string trim(const std::string &s) {
    size_t b = 0;
    while (b < s.size() && std::isspace((unsigned char)s[b])) {
        ++b;
    }

    size_t e = s.size();
    while (e > b && std::isspace((unsigned char)s[e - 1])) {
        --e;
    }

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
        } //endfor

        return out;
    } //endif

    return s;
}

static int find_col_optional(const std::vector<std::string> &header,
                             const std::string &target) {
    for (size_t i = 0; i < header.size(); ++i) {
        if (trim(unquote(header[i])) == target) {
            return (int)i;
        } //endif
    } //endfor

    return -1;
}

static int require_col(const std::vector<std::string> &header,
                       const std::string &target) {
    int idx = find_col_optional(header, target);
    if (idx < 0) {
        throw std::runtime_error("Missing required column: \"" + target + "\"");
    } //endif

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
        } //endif
    } //endwhile

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
    } //endif

    if (s.front() == '(' && s.back() == ')') {
        s = s.substr(1, s.size() - 2);
        s = trim(s);

        if (s.empty()) {
            return out;
        } //endif
    } //endif

    std::vector<std::string> parts;
    std::string token;

    for (size_t i = 0; i < s.size(); ++i) {
        const char c = s[i];

        if (c == ',') {
            parts.push_back(trim(token));
            token.clear();
        } else {
            token.push_back(c);
        } //endif
    } //endfor

    parts.push_back(trim(token));

    auto to_double_or_zero = [](const std::string &str) -> double {
        if (str.empty()) {
            return 0.0;
        } //endif

        return std::atof(str.c_str());
    };

    if (!parts.empty()) {
        out.value = to_double_or_zero(parts[0]);
    } //endif

    if (parts.size() > 1U) {
        out.stat = to_double_or_zero(parts[1]);
    } //endif

    if (parts.size() > 2U) {
        out.sys = to_double_or_zero(parts[2]);
    } //endif

    return out;
}

static bool finite_pos(double x) {
    return std::isfinite(x) && x > 0.0;
}

static Helicity helicity_from_string(const std::string &h) {
    if (h == "pos") {
        return Helicity::Plus;
    } //endif

    if (h == "neg") {
        return Helicity::Minus;
    } //endif

    return Helicity::Unpol;
}

static double beam_energy_for_label(const std::string &label) {
    if (label == "Sp19 Inb" || label == "10.2 GeV") {
        return 10.2;
    } //endif

    return 10.6;
}

static double wrap_phi_0_360(double phi_deg) {
    if (!std::isfinite(phi_deg)) {
        return std::numeric_limits<double>::quiet_NaN();
    } //endif

    double phi = std::fmod(phi_deg, 360.0);

    if (phi < 0.0) {
        phi += 360.0;
    } //endif

    return phi;
}

static double min_dist_to_0_or_360(double phi_deg) {
    const double phi = wrap_phi_0_360(phi_deg);

    if (!std::isfinite(phi)) {
        return std::numeric_limits<double>::infinity();
    } //endif

    return std::min(phi, 360.0 - phi);
}

static bool passes_edge_window(double dist) {
    if (!std::isfinite(dist)) {
        return false;
    } //endif

    if (kRequirePositiveDEdge && !(dist > 0.0)) {
        return false;
    } //endif

    if (dist > kMaxDEdgeForNormalizationDeg) {
        return false;
    } //endif

    return true;
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
            // Drop other punctuation deterministically.
        } //endif
    } //endfor

    std::string collapsed;
    collapsed.reserve(out.size());

    bool prev_us = false;
    for (size_t i = 0; i < out.size(); ++i) {
        if (out[i] == '_') {
            if (!prev_us) {
                collapsed.push_back('_');
            } //endif

            prev_us = true;
        } else {
            collapsed.push_back(out[i]);
            prev_us = false;
        } //endif
    } //endfor

    if (!collapsed.empty() && collapsed.back() == '_') {
        collapsed.pop_back();
    } //endif

    return collapsed;
}

static void ensure_output_dir_or_throw(const std::string &dir) {
    std::error_code ec;

    if (std::filesystem::exists(dir, ec)) {
        if (ec) {
            throw std::runtime_error("Failed to stat output directory: " + dir);
        } //endif

        return;
    } //endif

    if (!std::filesystem::create_directories(dir, ec)) {
        if (ec) {
            throw std::runtime_error("Failed to create output directory: " + dir);
        } //endif
    } //endif
}

// -----------------------------------------------------------------------------
// CSV write helpers
// -----------------------------------------------------------------------------

static std::string csv_join_fields(const std::vector<std::string> &fields) {
    std::ostringstream oss;

    for (size_t i = 0; i < fields.size(); ++i) {
        if (i) {
            oss << ",";
        } //endif

        oss << fields[i];
    } //endfor

    return oss.str();
}

static void write_lines_atomic_or_throw(const std::string &path,
                                        const std::vector<std::string> &lines) {
    std::filesystem::path p(path);
    std::filesystem::path dir = p.parent_path();

    std::string tmp = (dir / (p.filename().string() + ".tmp")).string();

    {
        std::ofstream ofs(tmp.c_str(), std::ios::out | std::ios::trunc);

        if (!ofs) {
            throw std::runtime_error("Failed to open temp CSV for write: " + tmp);
        } //endif

        for (size_t i = 0; i < lines.size(); ++i) {
            ofs << lines[i] << "\n";
        } //endfor

        ofs.close();

        if (!ofs) {
            throw std::runtime_error("Failed while writing temp CSV: " + tmp);
        } //endif
    }

    std::error_code ec;
    std::filesystem::rename(tmp, path, ec);

    if (ec) {
        std::filesystem::remove(path, ec);
        ec.clear();

        std::filesystem::rename(tmp, path, ec);

        if (ec) {
            throw std::runtime_error("Failed to rename temp CSV into place: " + tmp + " -> " + path);
        } //endif
    } //endif
}

static bool try_parse_double_blank_is_missing(const std::string &cell, double &out) {
    std::string s = trim(unquote(cell));

    if (s.empty()) {
        return false;
    } //endif

    out = std::atof(s.c_str());
    return std::isfinite(out);
}

static std::string format_double_for_csv(double x) {
    std::ostringstream oss;

    oss.setf(std::ios::fixed);
    oss << std::setprecision(6) << x;

    return oss.str();
}

// -----------------------------------------------------------------------------
// BH/model diagnostic grouping
// -----------------------------------------------------------------------------

enum class BhKmGroup {
    G_95_105 = 0,
    G_90_95_OR_105_110,
    G_85_90_OR_110_115,
    G_OUTSIDE_15PCT,
    G_INVALID
};

static BhKmGroup categorize_bh_over_model(double r) {
    if (!std::isfinite(r) || r <= 0.0) {
        return BhKmGroup::G_INVALID;
    } //endif

    if (r >= 0.95 && r <= 1.05) {
        return BhKmGroup::G_95_105;
    } //endif

    if ((r >= 0.90 && r < 0.95) || (r > 1.05 && r <= 1.10)) {
        return BhKmGroup::G_90_95_OR_105_110;
    } //endif

    if ((r > 0.85 && r < 0.90) || (r > 1.10 && r < 1.15)) {
        return BhKmGroup::G_85_90_OR_110_115;
    } //endif

    if (r <= 0.85 || r >= 1.15) {
        return BhKmGroup::G_OUTSIDE_15PCT;
    } //endif

    return BhKmGroup::G_INVALID;
}

struct GroupStyle {
    BhKmGroup group;
    int marker_style;
    int marker_color;
    std::string label;
};

static std::string ratio_name_for_metric(BhGroupMetric m) {
    if (m == BhGroupMetric::VGG) {
        return "BH/VGG";
    } //endif

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
        if (styles[i].group == g) {
            return (int)i;
        } //endif
    } //endfor

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
// Compute linear fit used for norm write
//
// This intentionally uses all points passed in. BH/VGG and BH/KM15 are not used
// to select points for this fit.
// -----------------------------------------------------------------------------

static bool compute_linear_fit_pol1_all_points(const std::vector<DepPoint> &pts,
                                               double x_min,
                                               double x_max,
                                               double &p0,
                                               double &p1,
                                               double &e0,
                                               double &e1,
                                               int &n_used) {
    p0 = 0.0;
    p1 = 0.0;
    e0 = 0.0;
    e1 = 0.0;
    n_used = 0;

    TGraphErrors gr;
    gr.SetName("gr_norm_fit_internal_all_edge_points");

    int n = 0;

    for (size_t i = 0; i < pts.size(); ++i) {
        if (!std::isfinite(pts[i].x) || !std::isfinite(pts[i].y)) {
            continue;
        } //endif

        const double ey = (std::isfinite(pts[i].ey) && pts[i].ey >= 0.0) ? pts[i].ey : 0.0;

        gr.SetPoint(n, pts[i].x, pts[i].y);
        gr.SetPointError(n, 0.0, ey);
        n += 1;
    } //endfor

    n_used = gr.GetN();

    if (n_used < 2) {
        return false;
    } //endif

    TF1 f("f_norm_pol1_internal_all_edge_points", "pol1", x_min, x_max);
    gr.Fit(&f, "Q0");

    p0 = f.GetParameter(0);
    p1 = f.GetParameter(1);
    e0 = f.GetParError(0);
    e1 = f.GetParError(1);

    if (!std::isfinite(p0) || !std::isfinite(p1)) {
        return false;
    } //endif

    return true;
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
    pad->SetLogx(0);
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
    tl.DrawLatex(0.14, 0.93, Form("%s, %s   (%s diagnostic)", label.c_str(), helicity.c_str(), ratio_name.c_str()));

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

        if (si < 0) {
            continue;
        } //endif

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
        } //endif
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
        } //endif
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

    double x_min = 0.0;
    double x_max = kMaxDEdgeForNormalizationDeg * 1.05;

    if (x_max <= x_min) {
        x_max = 10.0;
    } //endif

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
    tl.DrawLatex(0.12, 0.93, Form("%s, %s   (%s diagnostic)", label.c_str(), helicity.c_str(), ratio_name.c_str()));

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

        if (si < 0) {
            continue;
        } //endif

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
        } //endif
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
        } //endif
    } //endfor

    leg->Draw();
    leg->SetBit(kCanDelete);

    p->Update();
}

static void draw_all_edge_points_fit_dep_pad(TPad *p,
                                             const std::vector<DepPoint> &pts,
                                             const std::string &label,
                                             const std::string &helicity,
                                             const std::string &x_title,
                                             double x_min,
                                             double x_max,
                                             const std::string &fit_label) {
    p->SetLeftMargin(0.13);
    p->SetRightMargin(0.05);
    p->SetBottomMargin(0.14);
    p->SetTopMargin(0.12);
    p->SetGrid(1, 1);

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
    tl.DrawLatex(0.12, 0.93, Form("%s, %s   all edge points", label.c_str(), helicity.c_str()));

    TGraphErrors *gr = new TGraphErrors();
    gr->SetName(Form("gr_all_edge_fit_%s", sanitize_for_filename(fit_label).c_str()));
    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(kBlack);
    gr->SetLineColor(kBlack);
    gr->SetLineWidth(1);
    gr->SetMarkerSize(0.85);

    int n = 0;

    for (size_t i = 0; i < pts.size(); ++i) {
        if (!std::isfinite(pts[i].x) || !std::isfinite(pts[i].y)) {
            continue;
        } //endif

        gr->SetPoint(n, pts[i].x, pts[i].y);
        gr->SetPointError(n, 0.0, pts[i].ey);
        n += 1;
    } //endfor

    if (gr->GetN() <= 0) {
        delete gr;

        TLatex t0;
        t0.SetNDC(kTRUE);
        t0.SetTextSize(0.050);
        t0.DrawLatex(0.18, 0.70, "No points");

        p->Update();
        return;
    } //endif

    gr->Draw("PE1 SAME");
    gr->SetBit(kCanDelete);

    if (gr->GetN() >= 2) {
        TF1 *f_lin = new TF1(Form("f_lin_all_edge_%s", sanitize_for_filename(fit_label).c_str()), "pol1", x_min, x_max);

        gr->Fit(f_lin, "Q0");

        f_lin->SetLineColor(kRed+1);
        f_lin->SetLineStyle(2);
        f_lin->SetLineWidth(2);

        f_lin->Draw("SAME");
        f_lin->SetBit(kCanDelete);

        const double p0 = f_lin->GetParameter(0);
        const double p1 = f_lin->GetParameter(1);
        const double e0 = f_lin->GetParError(0);
        const double e1 = f_lin->GetParError(1);

        TLegend *legfit = new TLegend(0.54, 0.70, 0.93, 0.88);
        legfit->SetBorderSize(1);
        legfit->SetFillStyle(1001);
        legfit->SetFillColor(kWhite);
        legfit->SetTextSize(0.027);

        legfit->AddEntry((TObject*)0, "Linear fit (pol1)", "");
        legfit->AddEntry((TObject*)0, Form("p0 = %.4f #pm %.4f", p0, e0), "");
        legfit->AddEntry((TObject*)0, Form("p1 = %.4f #pm %.4f", p1, e1), "");
        legfit->AddEntry((TObject*)0, Form("N = %d", gr->GetN()), "");

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
        draw_all_edge_points_fit_dep_pad(p, pts, label, helicity, x_title, x_min, x_max, file_tag);
    }

    {
        TPad *p = (TPad*)c->cd(3);
        draw_all_groups_dep_pad(p, pts, label, helicity, x_title, x_min, x_max, BhGroupMetric::VGG);
    }

    {
        TPad *p = (TPad*)c->cd(4);
        draw_all_edge_points_fit_dep_pad(p, pts, label, helicity, x_title, x_min, x_max, file_tag);
    }

    const std::string out_png = out_dir + "/norm_xs_over_bh_vs_" + file_tag + "_" + label_tag + "_" + hel_tag + ".png";
    c->SaveAs(out_png.c_str());

    delete c;
}

// -----------------------------------------------------------------------------
// Row container
// -----------------------------------------------------------------------------

struct RowData {
    double xb_c;
    double q2_c;
    double t_c;
    double phi_deg;
    double dist_to_edge;
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

} // end anonymous namespace

static bool run_bh_normalization_study_impl(const std::string &csv_path,
                                            const std::string &label,
                                            const std::string &helicity) {
    try {
        std::ifstream ifs(csv_path.c_str());

        if (!ifs) {
            std::cerr << "[overall_norm] ERROR: cannot open CSV: " << csv_path << "\n";
            return false;
        } //endif

        std::vector<std::string> lines;
        std::string line;

        while (std::getline(ifs, line)) {
            lines.push_back(line);
        } //endwhile

        ifs.close();

        if (lines.empty()) {
            std::cerr << "[overall_norm] ERROR: CSV is empty: " << csv_path << "\n";
            return false;
        } //endif

        std::vector<std::string> header = split_csv_line(lines[0]);

        const std::string col_norm = "norm, " + label;
        const int c_norm = find_col_optional(header, col_norm);

        if (!normalization_helicity_allowed(label, helicity)) {
            std::cout << "[overall_norm] Skipping normalization study for unavailable helicity column: "
                      << label << ", " << helicity << "\n";
            return true;
        }

        // ---------------------------------------------------------------------
        // Fast bypass mode: write unity normalization and return immediately.
        // This intentionally does not require any other columns to exist.
        // ---------------------------------------------------------------------
        if (kSkipAllWorkWriteUnityNorm) {
            std::cout << "============================================================\n";
            std::cout << "[overall_norm] BH normalization study (UNITY OVERRIDE)\n";
            std::cout << "[overall_norm] CSV      : " << csv_path << "\n";
            std::cout << "[overall_norm] label    : " << label << "\n";
            std::cout << "[overall_norm] helicity : " << helicity << "\n";
            std::cout << "[overall_norm] mode     : SKIP ALL WORK; write norm=1.00 everywhere\n";

            if (c_norm < 0) {
                std::cout << "[overall_norm] NOTE: Column \"" << col_norm << "\" not found; nothing to write.\n";
                std::cout << "============================================================\n";
                return true;
            } //endif

            int n_rows_updated = 0;

            for (size_t r = 1; r < lines.size(); ++r) {
                if (lines[r].empty()) {
                    continue;
                } //endif

                std::vector<std::string> fields = split_csv_line(lines[r]);

                if (fields.size() != header.size()) {
                    continue;
                } //endif

                fields[(size_t)c_norm] = "1.00";
                lines[r] = csv_join_fields(fields);

                n_rows_updated += 1;
            } //endfor

            if (n_rows_updated > 0) {
                write_lines_atomic_or_throw(csv_path, lines);
                std::cout << "[overall_norm] Wrote unity norm values to CSV: " << n_rows_updated << " rows updated.\n";
            } else {
                std::cout << "[overall_norm] NOTE: no rows updated.\n";
            } //endif

            std::cout << "[overall_norm] Done (unity override).\n";
            std::cout << "============================================================\n";

            return true;
        } //endif

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

        const std::string col_e_theta = "e_theta, " + label;
        const std::string col_p_theta = "p_theta, " + label;
        const std::string col_g_theta = "g_theta, " + label;

        const int c_e_theta = require_col(header, col_e_theta);
        const int c_p_theta = require_col(header, col_p_theta);
        const int c_g_theta = require_col(header, col_g_theta);

        const std::string col_xs = "cross sections, ep->epg, exp, " + label + ", " + helicity;
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
                  << (kUseAllPointsWithinEdgeWindow ? "all points within edge window" : "closest-to-edge within edge window per kinematic cell")
                  << "\n";
        std::cout << "[overall_norm] edge cut : "
                  << (kRequirePositiveDEdge ? "0 < d_edge <= " : "d_edge <= ")
                  << std::fixed << std::setprecision(1) << kMaxDEdgeForNormalizationDeg
                  << " deg\n";
        std::cout << "[overall_norm] model ratios : diagnostic only; not used for point selection or norm fit\n";
        std::cout << "[overall_norm] norm x   : "
                  << (kNormXAxis == NormXAxis::PTheta ? "p_theta" : "xBavg")
                  << "\n";
        std::cout << "------------------------------------------------------------\n";

        std::vector<RowData> selected_rows;
        selected_rows.reserve(lines.size() > 1 ? lines.size() - 1 : 0);

        if (kUseAllPointsWithinEdgeWindow) {
            int n_candidate_rows = 0;
            int n_failed_kinematics = 0;
            int n_failed_xs = 0;
            int n_failed_edge = 0;

            for (size_t r = 1; r < lines.size(); ++r) {
                if (lines[r].empty()) {
                    continue;
                } //endif

                std::vector<std::string> fields = split_csv_line(lines[r]);

                if (fields.size() != header.size()) {
                    continue;
                } //endif

                n_candidate_rows += 1;

                const double xb_c = std::atof(trim(unquote(fields[(size_t)c_xbavg])).c_str());
                const double q2_c = std::atof(trim(unquote(fields[(size_t)c_q2avg])).c_str());
                const double t_c  = std::atof(trim(unquote(fields[(size_t)c_tavg])).c_str());
                const double phi  = std::atof(trim(unquote(fields[(size_t)c_phiavg])).c_str());

                const double e_th = std::atof(trim(unquote(fields[(size_t)c_e_theta])).c_str());
                const double p_th = std::atof(trim(unquote(fields[(size_t)c_p_theta])).c_str());
                const double g_th = std::atof(trim(unquote(fields[(size_t)c_g_theta])).c_str());

                if (!finite_pos(xb_c) || !finite_pos(q2_c) || !finite_pos(t_c) || !std::isfinite(phi)) {
                    n_failed_kinematics += 1;
                    continue;
                } //endif

                TripleCell xs = parse_tuple3(fields[(size_t)c_xs]);

                if (!(xs.value > 0.0) || !std::isfinite(xs.value)) {
                    n_failed_xs += 1;
                    continue;
                } //endif

                const double dist = min_dist_to_0_or_360(phi);

                if (!passes_edge_window(dist)) {
                    n_failed_edge += 1;
                    continue;
                } //endif

                RowData rd;
                rd.xb_c = xb_c;
                rd.q2_c = q2_c;
                rd.t_c  = t_c;
                rd.phi_deg = wrap_phi_0_360(phi);
                rd.dist_to_edge = dist;
                rd.xs = xs.value;
                rd.xs_stat = xs.stat;
                rd.e_theta = e_th;
                rd.p_theta = p_th;
                rd.g_theta = g_th;
                rd.csv_row_index = r;

                selected_rows.push_back(rd);
            } //endfor

            std::cout << "[overall_norm] Candidate CSV rows with matching columns : " << n_candidate_rows << "\n";
            std::cout << "[overall_norm] Rejected invalid kinematics             : " << n_failed_kinematics << "\n";
            std::cout << "[overall_norm] Rejected invalid/nonpositive xs         : " << n_failed_xs << "\n";
            std::cout << "[overall_norm] Rejected outside edge window            : " << n_failed_edge << "\n";

            if (selected_rows.empty()) {
                std::cerr << "[overall_norm] WARNING: no selected bins found within edge window.\n";
                return true;
            } //endif

            std::cout << "[overall_norm] Selected rows: " << selected_rows.size()
                      << " (all valid points within edge window)\n\n";
        } else {
            std::map<std::string, RowData> best;

            int n_candidate_rows = 0;
            int n_failed_kinematics = 0;
            int n_failed_xs = 0;
            int n_failed_edge = 0;

            for (size_t r = 1; r < lines.size(); ++r) {
                if (lines[r].empty()) {
                    continue;
                } //endif

                std::vector<std::string> fields = split_csv_line(lines[r]);

                if (fields.size() != header.size()) {
                    continue;
                } //endif

                n_candidate_rows += 1;

                const std::string xbmin_s = trim(unquote(fields[(size_t)c_xbmin]));
                const std::string xbmax_s = trim(unquote(fields[(size_t)c_xbmax]));
                const std::string q2min_s = trim(unquote(fields[(size_t)c_q2min]));
                const std::string q2max_s = trim(unquote(fields[(size_t)c_q2max]));
                const std::string tmin_s  = trim(unquote(fields[(size_t)c_tmin]));
                const std::string tmax_s  = trim(unquote(fields[(size_t)c_tmax]));

                const std::string key = cell_key_for_kin_bin(xbmin_s, xbmax_s, q2min_s, q2max_s, tmin_s, tmax_s);

                const double xb_c = std::atof(trim(unquote(fields[(size_t)c_xbavg])).c_str());
                const double q2_c = std::atof(trim(unquote(fields[(size_t)c_q2avg])).c_str());
                const double t_c  = std::atof(trim(unquote(fields[(size_t)c_tavg])).c_str());
                const double phi  = std::atof(trim(unquote(fields[(size_t)c_phiavg])).c_str());

                const double e_th = std::atof(trim(unquote(fields[(size_t)c_e_theta])).c_str());
                const double p_th = std::atof(trim(unquote(fields[(size_t)c_p_theta])).c_str());
                const double g_th = std::atof(trim(unquote(fields[(size_t)c_g_theta])).c_str());

                if (!finite_pos(xb_c) || !finite_pos(q2_c) || !finite_pos(t_c) || !std::isfinite(phi)) {
                    n_failed_kinematics += 1;
                    continue;
                } //endif

                TripleCell xs = parse_tuple3(fields[(size_t)c_xs]);

                if (!(xs.value > 0.0) || !std::isfinite(xs.value)) {
                    n_failed_xs += 1;
                    continue;
                } //endif

                const double dist = min_dist_to_0_or_360(phi);

                if (!passes_edge_window(dist)) {
                    n_failed_edge += 1;
                    continue;
                } //endif

                std::map<std::string, RowData>::iterator it = best.find(key);

                if (it == best.end() || dist < it->second.dist_to_edge) {
                    RowData rd;
                    rd.xb_c = xb_c;
                    rd.q2_c = q2_c;
                    rd.t_c  = t_c;
                    rd.phi_deg = wrap_phi_0_360(phi);
                    rd.dist_to_edge = dist;
                    rd.xs = xs.value;
                    rd.xs_stat = xs.stat;
                    rd.e_theta = e_th;
                    rd.p_theta = p_th;
                    rd.g_theta = g_th;
                    rd.csv_row_index = r;

                    best[key] = rd;
                } //endif
            } //endfor

            for (std::map<std::string, RowData>::const_iterator it = best.begin(); it != best.end(); ++it) {
                selected_rows.push_back(it->second);
            } //endfor

            std::cout << "[overall_norm] Candidate CSV rows with matching columns : " << n_candidate_rows << "\n";
            std::cout << "[overall_norm] Rejected invalid kinematics             : " << n_failed_kinematics << "\n";
            std::cout << "[overall_norm] Rejected invalid/nonpositive xs         : " << n_failed_xs << "\n";
            std::cout << "[overall_norm] Rejected outside edge window            : " << n_failed_edge << "\n";

            if (selected_rows.empty()) {
                std::cerr << "[overall_norm] WARNING: no selected bins found within edge window.\n";
                return true;
            } //endif

            std::cout << "[overall_norm] Selected rows: " << selected_rows.size()
                      << " (one per kinematic cell, closest to edge, within edge window)\n\n";
        } //endif

        std::cout
            << std::setw(8)  << "xB"
            << std::setw(10) << "Q2"
            << std::setw(12) << "t_abs"
            << std::setw(10) << "phi"
            << std::setw(10) << "d_edge"
            << std::setw(14) << "xs"
            << std::setw(14) << "xs/BH"
            << std::setw(14) << "BH/VGG"
            << std::setw(14) << "BH/KM15"
            << "\n";

        std::cout << std::string(8 + 10 + 12 + 10 + 10 + 14 * 4, '-') << "\n";

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

        double sumw_all = 0.0;
        double sumwx_all = 0.0;
        int n_weighted_used_all = 0;

        int n_bh_km15_95_105_total = 0;
        int n_bh_vgg_95_105_total = 0;

        for (size_t irow = 0; irow < selected_rows.size(); ++irow) {
            const RowData &br = selected_rows[irow];

            const double xb   = br.xb_c;
            const double q2   = br.q2_c;
            const double tpos = br.t_c;
            const double phi  = br.phi_deg;

            const double bh  = vgg_bh_only(xb, q2, tpos, phi, Ebeam);
            const double vgg = vgg_xs(xb, q2, tpos, phi, Ebeam, hel);
            const double km  = km15_xs(xb, q2, tpos, phi, Ebeam, hel);

            double xs_over_bh = 0.0;
            double xs_over_bh_err = 0.0;
            double bh_over_vgg = 0.0;
            double bh_over_km = 0.0;

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

            if (std::isfinite(bh_over_km) && bh_over_km >= 0.95 && bh_over_km <= 1.05) {
                n_bh_km15_95_105_total += 1;
            } //endif

            if (std::isfinite(bh_over_vgg) && bh_over_vgg >= 0.95 && bh_over_vgg <= 1.05) {
                n_bh_vgg_95_105_total += 1;
            } //endif

            std::cout
                << std::setw(8)  << std::fixed << std::setprecision(3) << xb
                << std::setw(10) << std::fixed << std::setprecision(2) << q2
                << std::setw(12) << std::fixed << std::setprecision(3) << tpos
                << std::setw(10) << std::fixed << std::setprecision(1) << phi
                << std::setw(10) << std::fixed << std::setprecision(1) << br.dist_to_edge
                << std::setw(14) << std::scientific << std::setprecision(3) << br.xs
                << std::setw(14) << std::fixed << std::setprecision(3) << xs_over_bh
                << std::setw(14) << std::fixed << std::setprecision(3) << bh_over_vgg
                << std::setw(14) << std::fixed << std::setprecision(3) << bh_over_km
                << "\n";

            if (std::isfinite(xs_over_bh) && xs_over_bh >= 0.0 &&
                std::isfinite(bh_over_km) && bh_over_km > 0.0 &&
                std::isfinite(bh_over_vgg) && bh_over_vgg > 0.0) {

                PlotPoint p;
                p.d_edge_deg      = br.dist_to_edge;
                p.xs_over_bh      = xs_over_bh;
                p.xs_over_bh_err  = xs_over_bh_err;
                p.bh_over_km15    = bh_over_km;
                p.bh_over_vgg     = bh_over_vgg;

                plot_pts_dedge.push_back(p);

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

                if (std::isfinite(xs_over_bh_err) && xs_over_bh_err > 0.0) {
                    const double w = 1.0 / (xs_over_bh_err * xs_over_bh_err);

                    sumw_all  += w;
                    sumwx_all += w * xs_over_bh;

                    n_weighted_used_all += 1;
                } //endif
            } //endif
        } //endfor

        const std::string out_root = kOutputRoot;
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
                                 "Q^{2} (GeV^{2})",
                                 "q2",
                                 1.0, 6.0,
                                 dep_q2);

        make_dependence_plot_2x2(out_dir, label, helicity,
                                 "-t (GeV^{2})",
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

        // ---------------------------------------------------------------------
        // Write "norm, <label>" by evaluating the fitted line vs chosen variable.
        //
        // The fit now uses all selected edge-window points. It does not require
        // BH/KM15 or BH/VGG to lie in [0.95, 1.05].
        //
        // If kNormXAxis == PTheta:
        //   - Fit uses dep_p_theta
        //   - Writes per-row using "p_theta, <label>"
        //
        // If kNormXAxis == XB:
        //   - Fit uses dep_xb
        //   - Writes per-row using "xBavg, <label>"
        // ---------------------------------------------------------------------

        if (c_norm < 0) {
            std::cout << "\n";
            std::cout << "[overall_norm] NOTE: Column \"" << col_norm << "\" not found; skipping norm write.\n";
        } else {
            const std::vector<DepPoint> *dep_for_fit = nullptr;
            int c_x_for_write = -1;
            double fit_xmin = 0.0;
            double fit_xmax = 0.0;
            std::string x_name;

            if (kNormXAxis == NormXAxis::PTheta) {
                dep_for_fit = &dep_p_theta;
                c_x_for_write = c_p_theta;
                fit_xmin = 0.0;
                fit_xmax = 80.0;
                x_name = "p_theta";
            } else {
                dep_for_fit = &dep_xb;
                c_x_for_write = c_xbavg;
                fit_xmin = 0.0;
                fit_xmax = 0.6;
                x_name = "xBavg";
            } //endif

            double fit_p0 = 0.0;
            double fit_p1 = 0.0;
            double fit_e0 = 0.0;
            double fit_e1 = 0.0;
            int n_fit_used = 0;

            const bool fit_ok = compute_linear_fit_pol1_all_points(*dep_for_fit,
                                                                    fit_xmin,
                                                                    fit_xmax,
                                                                    fit_p0,
                                                                    fit_p1,
                                                                    fit_e0,
                                                                    fit_e1,
                                                                    n_fit_used);

            std::cout << "\n";

            if (!fit_ok) {
                std::cout << "[overall_norm] WARNING: edge-window fit failed; skipping norm write.\n";
                std::cout << "              need >=2 valid edge-window points\n";
                std::cout << "              norm x = " << x_name << "\n";
            } else {
                std::cout << "[overall_norm] Writing \"" << col_norm << "\" using pol1 fit vs " << x_name << ":\n";
                std::cout << "              p0 = " << std::fixed << std::setprecision(6) << fit_p0
                          << " +/- " << std::fixed << std::setprecision(6) << fit_e0 << "\n";
                std::cout << "              p1 = " << std::fixed << std::setprecision(6) << fit_p1
                          << " +/- " << std::fixed << std::setprecision(6) << fit_e1 << "\n";
                std::cout << "              edge-window points used = " << n_fit_used << "\n";
                std::cout << "              model ratios were diagnostic only\n";

                int n_rows_written = 0;
                int n_rows_skipped_blank_x = 0;

                for (size_t r = 1; r < lines.size(); ++r) {
                    if (lines[r].empty()) {
                        continue;
                    } //endif

                    std::vector<std::string> fields = split_csv_line(lines[r]);

                    if (fields.size() != header.size()) {
                        continue;
                    } //endif

                    double x_val = std::numeric_limits<double>::quiet_NaN();

                    if (!try_parse_double_blank_is_missing(fields[(size_t)c_x_for_write], x_val)) {
                        n_rows_skipped_blank_x += 1;
                        continue;
                    } //endif

                    const double norm_val = fit_p0 + fit_p1 * x_val;

                    if (!std::isfinite(norm_val)) {
                        continue;
                    } //endif

                    fields[(size_t)c_norm] = format_double_for_csv(norm_val);
                    lines[r] = csv_join_fields(fields);

                    n_rows_written += 1;
                } //endfor

                if (n_rows_written > 0) {
                    write_lines_atomic_or_throw(csv_path, lines);

                    std::cout << "[overall_norm] Wrote norm values to CSV: " << n_rows_written << " rows updated.\n";
                    std::cout << "[overall_norm] Skipped due to blank " << x_name << ": " << n_rows_skipped_blank_x << " rows.\n";
                } else {
                    std::cout << "[overall_norm] NOTE: no rows updated.\n";
                } //endif
            } //endif
        } //endif

        std::cout << "\n";
        std::cout << "------------------------------------------------------------\n";
        std::cout << "[overall_norm] Edge-window weighted xs/BH summary\n";
        std::cout << "[overall_norm] Selected rows                 : " << selected_rows.size() << "\n";
        std::cout << "[overall_norm] Valid plotted/fitted points    : " << dep_xb.size() << "\n";
        std::cout << "[overall_norm] Points with BH/KM15 in [0.95, 1.05] : " << n_bh_km15_95_105_total << "\n";
        std::cout << "[overall_norm] Points with BH/VGG  in [0.95, 1.05] : " << n_bh_vgg_95_105_total << "\n";
        std::cout << "[overall_norm] Points used in weighted mean (err>0): " << n_weighted_used_all << "\n";

        if (sumw_all > 0.0) {
            const double mean = sumwx_all / sumw_all;
            const double err  = std::sqrt(1.0 / sumw_all);

            std::cout << "[overall_norm] Weighted mean xs/BH over all edge points : "
                      << std::fixed << std::setprecision(6) << mean << "\n";
            std::cout << "[overall_norm] Weighted stat unc                       : "
                      << std::fixed << std::setprecision(6) << err << "\n";
        } else {
            std::cout << "[overall_norm] Weighted mean xs/BH over all edge points : N/A\n";
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

bool update_overall_normalization_study_csv(const std::string &csv_path,
                                            const std::string &label,
                                            const std::string &helicity,
                                            const OverallNormalizationOptions &options) {
    kSkipAllWorkWriteUnityNorm = options.override_to_unity;
    kUseAllPointsWithinEdgeWindow = options.use_all_points_within_edge_window;
    kRequirePositiveDEdge = options.require_positive_dedge;
    kMaxDEdgeForNormalizationDeg = options.max_dedge_for_normalization_deg;
    kNormXAxis = (options.norm_x_axis == OverallNormXAxis::PTheta) ? NormXAxis::PTheta : NormXAxis::XB;
    kOutputRoot = options.output_dir;

    return run_bh_normalization_study_impl(csv_path, label, helicity);
}

bool print_bh_normalization_study(const std::string &csv_path,
                                  const std::string &label,
                                  const std::string &helicity) {
    OverallNormalizationOptions options;
    options.override_to_unity = true;
    return update_overall_normalization_study_csv(csv_path, label, helicity, options);
}
