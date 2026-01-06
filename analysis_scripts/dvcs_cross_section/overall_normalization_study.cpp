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

// ROOT (legacy + new 1D dependence plots; 2D plot removed)
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TColor.h>
#include <TROOT.h>

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
// Grouping for the xs/BH plots: by BH/KM15 range (kept exactly as-is)
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

struct DepPoint {
    double x;
    double y;
    double ey;
    double bh_over_km15;
};

static void make_normalization_plots_legacy_1d_only(const std::string &out_dir,
                                                    const std::string &label,
                                                    const std::string &helicity,
                                                    const std::vector<PlotPoint> &pts) {
    ensure_output_dir_or_throw(out_dir);

    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);

    const std::string label_tag = sanitize_for_filename(label);
    const std::string hel_tag   = sanitize_for_filename(helicity);

    // x-axis (d_edge) starts at 1 degree
    const double x_min = 1.0;
    double x_max = x_min * 10.0;
    for (size_t i = 0; i < pts.size(); ++i) {
        if (std::isfinite(pts[i].d_edge_deg) && pts[i].d_edge_deg > x_max) {
            x_max = pts[i].d_edge_deg;
        }
    }
    x_max *= 1.10;

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
        std::cerr << "[overall_norm] WARNING: no points available for 1D plot.\n";
    }

    for (size_t si = 0; si < graphs.size(); ++si) {
        delete graphs[si];
    }
    delete c;
}

static void make_dependence_plot(const std::string &out_dir,
                                 const std::string &label,
                                 const std::string &helicity,
                                 const std::string &x_title,
                                 const std::string &plot_tag,
                                 const std::vector<DepPoint> &pts) {
    ensure_output_dir_or_throw(out_dir);

    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);

    const std::string label_tag = sanitize_for_filename(label);
    const std::string hel_tag   = sanitize_for_filename(helicity);

    const std::vector<GroupStyle> styles = make_group_styles();

    std::vector<TGraphErrors*> graphs(styles.size(), (TGraphErrors*)nullptr);
    std::vector<int> n_in(styles.size(), 0);

    double x_min = std::numeric_limits<double>::infinity();
    double x_max = -std::numeric_limits<double>::infinity();
    bool any_x = false;

    for (size_t si = 0; si < styles.size(); ++si) {
        graphs[si] = new TGraphErrors();
        graphs[si]->SetName(Form("gr_dep_%zu_%s", si, plot_tag.c_str()));
        graphs[si]->SetMarkerStyle(styles[si].marker_style);
        graphs[si]->SetMarkerColor(styles[si].marker_color);
        graphs[si]->SetLineColor(styles[si].marker_color);
        graphs[si]->SetLineWidth(1);
    }

    for (size_t i = 0; i < pts.size(); ++i) {
        if (!std::isfinite(pts[i].x)) continue;
        if (!std::isfinite(pts[i].y) || pts[i].y < 0.0) continue;

        const BhKmGroup g = categorize_bh_over_km15(pts[i].bh_over_km15);
        const int si = style_index_for_group(g, styles);
        if (si < 0) continue;

        const int n = n_in[(size_t)si];
        graphs[(size_t)si]->SetPoint(n, pts[i].x, pts[i].y);
        graphs[(size_t)si]->SetPointError(n, 0.0, (std::isfinite(pts[i].ey) && pts[i].ey >= 0.0 ? pts[i].ey : 0.0));
        n_in[(size_t)si] += 1;

        if (pts[i].x < x_min) x_min = pts[i].x;
        if (pts[i].x > x_max) x_max = pts[i].x;
        any_x = true;
    }

    if (!any_x || !(x_max > x_min)) {
        std::cerr << "[overall_norm] WARNING: no usable points for dependence plot: " << plot_tag << "\n";
        for (size_t si = 0; si < graphs.size(); ++si) {
            delete graphs[si];
        }
        return;
    }

    // Pad x-range a bit
    const double dx = x_max - x_min;
    x_min -= 0.05 * dx;
    x_max += 0.05 * dx;

    TCanvas *c = new TCanvas(Form("c_dep_%s", plot_tag.c_str()),
                             Form("c_dep_%s", plot_tag.c_str()),
                             950, 700);
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
        graphs[(size_t)first_nonempty]->GetXaxis()->SetTitle(x_title.c_str());
        graphs[(size_t)first_nonempty]->GetYaxis()->SetTitle("xs/BH");
        graphs[(size_t)first_nonempty]->GetXaxis()->SetLimits(x_min, x_max);
        graphs[(size_t)first_nonempty]->GetYaxis()->SetRangeUser(0.0, 3.0);

        TLatex tl;
        tl.SetNDC(kTRUE);
        tl.SetTextSize(0.040);
        tl.DrawLatex(0.14, 0.93, Form("Normalization study: xs/BH vs %s (%s, %s)",
                                      x_title.c_str(), label.c_str(), helicity.c_str()));

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
            out_dir + "/norm_1d_xs_over_bh_vs_" + plot_tag + "_" + label_tag + "_" + hel_tag + ".png";
        c->SaveAs(out_png.c_str());
    } else {
        std::cerr << "[overall_norm] WARNING: no points available for dependence plot: " << plot_tag << "\n";
    }

    for (size_t si = 0; si < graphs.size(); ++si) {
        delete graphs[si];
    }
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

        // Use per-period center columns
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

        // Best phi (closest to 0/360) per (xB,Q2,t) bin
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

            const std::string key_legacy = cell_key_for_kin_bin(xbmin_s, xbmax_s, q2min_s, q2max_s, tmin_s, tmax_s);

            const double xb_c   = std::atof(trim(unquote(fields[c_xbavg])).c_str());
            const double q2_c   = std::atof(trim(unquote(fields[c_q2avg])).c_str());
            const double t_c    = std::atof(trim(unquote(fields[c_tavg])).c_str()); // already positive t_abs
            const double phi    = std::atof(trim(unquote(fields[c_phiavg])).c_str());

            if (!finite_pos(xb_c) || !finite_pos(q2_c) || !finite_pos(t_c) || !std::isfinite(phi)) {
                continue;
            }

            TripleCell xs = parse_tuple3(fields[c_xs]);
            if (!(xs.value > 0.0) || !std::isfinite(xs.value)) {
                continue;
            }

            const double dist = min_dist_to_0_or_360(phi);
            std::map<std::string, BestPhiRow>::iterator it = best.find(key_legacy);
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
                best[key_legacy] = br;
            }
        }

        if (best.empty()) {
            std::cerr << "[overall_norm] WARNING: no bins found with positive cross section values.\n";
            return true;
        }

        std::cout << "[overall_norm] Unique kinematic bins found (legacy best-phi): " << best.size()
                  << " (from " << n_rows << " CSV data rows)\n\n";

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

        std::vector<PlotPoint> plot_pts;
        plot_pts.reserve(best.size());

        std::vector<DepPoint> dep_xb;
        std::vector<DepPoint> dep_q2;
        std::vector<DepPoint> dep_t;

        dep_xb.reserve(best.size());
        dep_q2.reserve(best.size());
        dep_t.reserve(best.size());

        double sumw = 0.0;
        double sumwx = 0.0;
        int n_weighted_used = 0;
        int n_in_95_105_total = 0;

        for (std::map<std::string, BestPhiRow>::const_iterator it = best.begin();
             it != best.end(); ++it) {

            const BestPhiRow &br = it->second;

            const double xb   = br.xb_c;
            const double q2   = br.q2_c;
            const double tpos = br.t_c;       // t_abs already positive
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

            if (std::isfinite(xs_over_bh) && xs_over_bh >= 0.0 &&
                std::isfinite(xs_over_bh_err) && xs_over_bh_err >= 0.0 &&
                std::isfinite(bh_over_km) && bh_over_km > 0.0) {

                if (std::isfinite(xb)) {
                    DepPoint d;
                    d.x = xb;
                    d.y = xs_over_bh;
                    d.ey = xs_over_bh_err;
                    d.bh_over_km15 = bh_over_km;
                    dep_xb.push_back(d);
                }

                if (std::isfinite(q2)) {
                    DepPoint d;
                    d.x = q2;
                    d.y = xs_over_bh;
                    d.ey = xs_over_bh_err;
                    d.bh_over_km15 = bh_over_km;
                    dep_q2.push_back(d);
                }

                if (std::isfinite(tpos)) {
                    DepPoint d;
                    d.x = tpos; // DO NOT NEGATE (t_abs already positive)
                    d.y = xs_over_bh;
                    d.ey = xs_over_bh_err;
                    d.bh_over_km15 = bh_over_km;
                    dep_t.push_back(d);
                }
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
        ensure_output_dir_or_throw(out_dir);

        make_normalization_plots_legacy_1d_only(out_dir, label, helicity, plot_pts);

        make_dependence_plot(out_dir, label, helicity, "x_{B}", "xb", dep_xb);
        make_dependence_plot(out_dir, label, helicity, "Q^{2} (GeV^{2})", "q2", dep_q2);
        make_dependence_plot(out_dir, label, helicity, "|t| (GeV^{2})", "t", dep_t);

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