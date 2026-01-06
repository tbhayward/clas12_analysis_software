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
#include <TH1F.h>
#include <TF1.h>
#include <TFitResultPtr.h>
#include <TFitResult.h>
#include <TMatrixDSym.h>
#include <TMath.h>

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
// Legacy grouping for the 1D xs/BH vs d_edge plot (kept exactly as-is here)
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

static void make_normalization_plots_legacy(const std::string &out_dir,
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

            if (g2->GetXaxis()) g2->GetXaxis()->SetLimits(0.1, x_max);

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

// -----------------------------------------------------------------------------
// New grid machinery: fit xs(phi) to A(1 + B cos + C cos2 + D cos3) and evaluate
// at phi0=0.001 deg.
// -----------------------------------------------------------------------------

struct PhiPoint {
    double phi_deg;
    double xs;
    double xs_stat;
};

struct KinKey {
    // Use the exact CSV edge strings for deterministic grouping
    std::string xbmin_s, xbmax_s;
    std::string q2min_s, q2max_s;
    std::string tmin_s,  tmax_s;

    bool operator<(const KinKey &o) const {
        if (xbmin_s != o.xbmin_s) return xbmin_s < o.xbmin_s;
        if (xbmax_s != o.xbmax_s) return xbmax_s < o.xbmax_s;
        if (q2min_s != o.q2min_s) return q2min_s < o.q2min_s;
        if (q2max_s != o.q2max_s) return q2max_s < o.q2max_s;
        if (tmin_s  != o.tmin_s)  return tmin_s  < o.tmin_s;
        return tmax_s < o.tmax_s;
    }
};

struct XbBinKey {
    std::string xbmin_s, xbmax_s;

    bool operator<(const XbBinKey &o) const {
        if (xbmin_s != o.xbmin_s) return xbmin_s < o.xbmin_s;
        return xbmax_s < o.xbmax_s;
    }
};

struct Q2BinKey {
    std::string q2min_s, q2max_s;

    bool operator<(const Q2BinKey &o) const {
        const double a0 = std::atof(q2min_s.c_str());
        const double b0 = std::atof(o.q2min_s.c_str());
        if (a0 != b0) return a0 < b0;
        const double a1 = std::atof(q2max_s.c_str());
        const double b1 = std::atof(o.q2max_s.c_str());
        if (a1 != b1) return a1 < b1;
        if (q2min_s != o.q2min_s) return q2min_s < o.q2min_s;
        return q2max_s < o.q2max_s;
    }
};

struct TBinKey {
    std::string tmin_s, tmax_s;

    bool operator<(const TBinKey &o) const {
        const double a0 = std::atof(tmin_s.c_str());
        const double b0 = std::atof(o.tmin_s.c_str());
        if (a0 != b0) return a0 < b0;
        const double a1 = std::atof(tmax_s.c_str());
        const double b1 = std::atof(o.tmax_s.c_str());
        if (a1 != b1) return a1 < b1;
        if (tmin_s != o.tmin_s) return tmin_s < o.tmin_s;
        return tmax_s < o.tmax_s;
    }
};

struct BinMeta {
    double xb_c;
    double q2_c;
    double t_c; // positive |t|
};

static std::string build_cos_series_formula(int max_harmonic) {
    // x is phi in degrees. Convert internally to radians.
    std::ostringstream oss;
    oss << "[0]*(1";
    for (int n = 1; n <= max_harmonic; ++n) {
        oss << " + [" << n << "]*cos(" << n << "*x*(TMath::Pi()/180.0))";
    }
    oss << ")";
    return oss.str();
}

static bool fit_cos_series(const std::vector<PhiPoint> &pts_in,
                           int &used_harmonic,
                           double &A, double &dA,
                           std::vector<double> &pars,
                           TMatrixDSym &cov_out,
                           std::string &status_msg) {
    status_msg = "";

    // Keep only points with finite phi and positive xs; for fitting weights require xs_stat>0.
    std::vector<PhiPoint> pts;
    pts.reserve(pts_in.size());
    for (size_t i = 0; i < pts_in.size(); ++i) {
        if (!std::isfinite(pts_in[i].phi_deg)) continue;
        if (!std::isfinite(pts_in[i].xs) || pts_in[i].xs <= 0.0) continue;
        if (!std::isfinite(pts_in[i].xs_stat) || pts_in[i].xs_stat <= 0.0) continue;
        pts.push_back(pts_in[i]);
    }

    if (pts.empty()) {
        status_msg = "No points with positive xs and xs_stat>0";
        return false;
    }

    // Sort by phi for nicer drawing later (fit doesn't care).
    std::sort(pts.begin(), pts.end(), [](const PhiPoint &a, const PhiPoint &b) {
        return a.phi_deg < b.phi_deg;
    });

    // Try harmonics 3 -> 2 -> 1 -> 0 until enough points.
    int H = 3;
    while (H >= 0) {
        const int npar = 1 + H;
        if ((int)pts.size() >= npar) break;
        --H;
    }

    if (H < 0) {
        status_msg = "Too few points to fit even A-only";
        return false;
    }

    // Build graph (stack object is fine here; not drawn/persisted)
    TGraphErrors gr;
    gr.SetName("gr_fit");
    gr.SetMarkerStyle(20);
    gr.SetMarkerColor(kBlack);
    gr.SetLineColor(kBlack);

    int n = 0;
    double ysum = 0.0;
    for (size_t i = 0; i < pts.size(); ++i) {
        gr.SetPoint(n, pts[i].phi_deg, pts[i].xs);
        gr.SetPointError(n, 0.0, pts[i].xs_stat);
        ysum += pts[i].xs;
        n += 1;
    }

    const double A_init = ysum / (double)pts.size();

    // Fit; if fit fails at H, reduce harmonic and retry.
    while (H >= 0) {
        const std::string form = build_cos_series_formula(H);
        TF1 f("f_cos", form.c_str(), 0.0, 360.0);
        f.SetNpx(720);

        // Parameter initialization
        f.SetParameter(0, A_init);
        for (int i = 1; i <= H; ++i) {
            f.SetParameter(i, 0.0);
        }

        // Quiet, return fit result, do not draw automatically
        TFitResultPtr r = gr.Fit(&f, "QSN");
        const int fit_status = (int)r;

        if (fit_status == 0 && r.Get() != nullptr && r->IsValid()) {
            used_harmonic = H;

            pars.clear();
            pars.resize((size_t)(1 + H), 0.0);
            for (int ip = 0; ip <= H; ++ip) {
                pars[(size_t)ip] = f.GetParameter(ip);
            }

            A  = pars[0];
            dA = f.GetParError(0);

            cov_out.ResizeTo(1 + H, 1 + H);
            cov_out = r->GetCovarianceMatrix();

            return true;
        }

        // Reduce harmonic and retry
        --H;
    }

    status_msg = "Fit failed for all harmonic orders";
    return false;
}

static bool eval_fit_and_unc_at_phi0(double phi0_deg,
                                     int used_harmonic,
                                     const std::vector<double> &pars,
                                     const TMatrixDSym &cov,
                                     double &xs_fit0,
                                     double &dxs_fit0) {
    xs_fit0 = std::numeric_limits<double>::quiet_NaN();
    dxs_fit0 = std::numeric_limits<double>::quiet_NaN();

    if (pars.empty()) return false;
    if (used_harmonic < 0) return false;
    if ((int)pars.size() != 1 + used_harmonic) return false;
    if (cov.GetNrows() != (int)pars.size()) return false;

    const double phi_rad = phi0_deg * (TMath::Pi() / 180.0);

    // f = A*(1 + sum_{n=1..H} p_n cos(n phi))
    double s = 1.0;
    for (int n = 1; n <= used_harmonic; ++n) {
        s += pars[(size_t)n] * std::cos((double)n * phi_rad);
    }
    xs_fit0 = pars[0] * s;

    // Gradient wrt parameters:
    // df/dA = (1 + sum p_n cos(n phi)) = s
    // df/dp_n = A*cos(n phi)
    std::vector<double> g;
    g.resize((size_t)(1 + used_harmonic), 0.0);
    g[0] = s;
    for (int n = 1; n <= used_harmonic; ++n) {
        g[(size_t)n] = pars[0] * std::cos((double)n * phi_rad);
    }

    double var = 0.0;
    for (int i = 0; i <= used_harmonic; ++i) {
        for (int j = 0; j <= used_harmonic; ++j) {
            var += g[(size_t)i] * cov(i, j) * g[(size_t)j];
        }
    }
    dxs_fit0 = (var > 0.0 ? std::sqrt(var) : 0.0);
    return std::isfinite(xs_fit0) && std::isfinite(dxs_fit0);
}

static std::string fmt_edge(const std::string &s) {
    // For display in TLatex: keep as-is but trim
    return trim(unquote(s));
}

static std::string fmt_edge_fname(const std::string &s) {
    // For filenames: sanitize deterministically
    return sanitize_for_filename(fmt_edge(s));
}

static std::string neg_edge_str(const std::string &s) {
    std::string a = fmt_edge(s);
    if (!a.empty() && a[0] == '-') return a;
    return "-" + a;
}

static void draw_pad_text_only(const std::string &line1,
                               const std::string &line2,
                               const std::string &line3) {
    TLatex tl;
    tl.SetNDC(kTRUE);
    tl.SetTextSize(0.060);
    tl.DrawLatex(0.12, 0.78, line1.c_str());
    tl.SetTextSize(0.055);
    tl.DrawLatex(0.12, 0.62, line2.c_str());
    if (!line3.empty()) {
        tl.DrawLatex(0.12, 0.46, line3.c_str());
    } //endif
}

// NOTE: Keep per-pad TGraphErrors / TF1 / TLegend objects alive until after SaveAs.
// ROOT pads store pointers to primitives; deleting them early results in blank pads.
static void make_normalization_grids(const std::string &out_dir_grid,
                                     const std::string &label,
                                     const std::string &helicity,
                                     double Ebeam,
                                     Helicity hel,
                                     const std::map<KinKey, std::vector<PhiPoint> > &bin_points,
                                     const std::map<KinKey, BinMeta> &bin_meta) {
    ensure_output_dir_or_throw(out_dir_grid);

    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);

    const double phi0_deg = 0.001;

    // Group bins by xB edges
    std::map<XbBinKey, std::set<Q2BinKey> > xb_to_q2;
    std::map<XbBinKey, std::set<TBinKey> >  xb_to_t;
    std::map<XbBinKey, std::set<KinKey> >   xb_to_bins;

    for (std::map<KinKey, std::vector<PhiPoint> >::const_iterator it = bin_points.begin();
         it != bin_points.end(); ++it) {

        const KinKey &k = it->first;
        XbBinKey xb{ k.xbmin_s, k.xbmax_s };
        Q2BinKey q2{ k.q2min_s, k.q2max_s };
        TBinKey  tt{ k.tmin_s,  k.tmax_s  };

        xb_to_q2[xb].insert(q2);
        xb_to_t[xb].insert(tt);
        xb_to_bins[xb].insert(k);
    }

    // Print header for per-bin fit summary
    std::cout << "------------------------------------------------------------\n";
    std::cout << "[overall_norm] Per-bin harmonic fits: xs(#phi) = A(1 + Bcos#phi + Ccos2#phi + Dcos3#phi)\n";
    std::cout << "[overall_norm] Evaluation at #phi = " << std::fixed << std::setprecision(3) << phi0_deg << " (deg)\n";
    std::cout << std::fixed << std::setprecision(5);
    std::cout
        << std::setw(10) << "xB"
        << std::setw(10) << "Q2"
        << std::setw(10) << "-t"
        << std::setw(12) << "A"
        << std::setw(12) << "dA"
        << std::setw(12) << "BH0"
        << std::setw(12) << "BH/KM15"
        << std::setw(12) << "xs_fit0"
        << std::setw(12) << "dxs_fit0"
        << std::setw(14) << "dxs/BH0"
        << "\n";
    std::cout << std::string(10+10+10+12+12+12+12+12+12+14, '-') << "\n";

    // Make one canvas per xB bin
    for (std::map<XbBinKey, std::set<KinKey> >::const_iterator it_xb = xb_to_bins.begin();
         it_xb != xb_to_bins.end(); ++it_xb) {

        const XbBinKey &xbk = it_xb->first;

        // Ordered lists of q2 and t for this xB bin
        std::vector<Q2BinKey> q2_bins;
        std::vector<TBinKey>  t_bins;

        {
            const std::set<Q2BinKey> &qs = xb_to_q2[xbk];
            q2_bins.assign(qs.begin(), qs.end());
        }
        {
            const std::set<TBinKey> &ts = xb_to_t[xbk];
            t_bins.assign(ts.begin(), ts.end());
        }

        const int nrows = (int)q2_bins.size();
        const int ncols = (int)t_bins.size();

        if (nrows <= 0 || ncols <= 0) {
            continue;
        } //endif

        // ---------------------------------------------------------------------
        // NEW: determine shared log-y range for this entire xB canvas
        // ---------------------------------------------------------------------
        double y_min_pos = std::numeric_limits<double>::infinity();
        double y_max_all = 0.0;
        bool any_points = false;

        for (int ir = 0; ir < nrows; ++ir) {
            for (int ic = 0; ic < ncols; ++ic) {
                const Q2BinKey &q2k = q2_bins[(size_t)ir];
                const TBinKey  &tk  = t_bins[(size_t)ic];

                KinKey kk;
                kk.xbmin_s = xbk.xbmin_s;
                kk.xbmax_s = xbk.xbmax_s;
                kk.q2min_s = q2k.q2min_s;
                kk.q2max_s = q2k.q2max_s;
                kk.tmin_s  = tk.tmin_s;
                kk.tmax_s  = tk.tmax_s;

                std::map<KinKey, std::vector<PhiPoint> >::const_iterator itp = bin_points.find(kk);
                if (itp == bin_points.end()) continue;

                const std::vector<PhiPoint> &pts = itp->second;
                for (size_t i = 0; i < pts.size(); ++i) {
                    if (!std::isfinite(pts[i].phi_deg)) continue;
                    if (!std::isfinite(pts[i].xs) || pts[i].xs <= 0.0) continue;

                    const double ey = (std::isfinite(pts[i].xs_stat) && pts[i].xs_stat > 0.0 ? pts[i].xs_stat : 0.0);

                    const double hi = pts[i].xs + ey;
                    if (hi > y_max_all) y_max_all = hi;

                    const double lo_candidate = pts[i].xs - ey;
                    const double lo = (lo_candidate > 0.0 ? lo_candidate : pts[i].xs);
                    if (lo > 0.0 && lo < y_min_pos) y_min_pos = lo;

                    any_points = true;
                } //endfor
            } //endfor
        } //endfor

        if (!any_points || !(y_max_all > 0.0) || !std::isfinite(y_min_pos) || !(y_min_pos > 0.0)) {
            std::cerr << "[overall_norm] WARNING: no usable positive xs points for xB canvas "
                      << "[" << fmt_edge(xbk.xbmin_s) << ", " << fmt_edge(xbk.xbmax_s) << "]. Skipping.\n";
            continue;
        } //endif

        // Pad the shared range; enforce strict positivity for log-scale
        double y_min_global = y_min_pos * 0.80;
        double y_max_global = y_max_all * 1.25;

        if (!(y_min_global > 0.0)) y_min_global = y_min_pos * 0.50;
        if (!(y_min_global > 0.0)) y_min_global = 1e-12;

        if (!(y_max_global > y_min_global)) y_max_global = y_min_global * 10.0;

        const int W = 320 * ncols + 220;
        const int H = 260 * nrows + 180;

        const std::string cname = "c_grid_" + fmt_edge_fname(xbk.xbmin_s) + "_" + fmt_edge_fname(xbk.xbmax_s);
        TCanvas *c = new TCanvas(cname.c_str(), cname.c_str(), W, H);

        // IMPORTANT: keep per-canvas objects alive through SaveAs
        std::vector<std::unique_ptr<TGraphErrors> > keep_graphs;
        std::vector<std::unique_ptr<TF1> >         keep_fits;
        std::vector<std::unique_ptr<TLegend> >     keep_legs;

        // Title pad
        TPad *ptitle = new TPad((cname + "_title").c_str(), "title", 0.0, 0.92, 1.0, 1.0);
        ptitle->SetFillStyle(0);
        ptitle->SetFrameFillStyle(0);
        ptitle->Draw();
        ptitle->cd();

        TLatex tlt;
        tlt.SetNDC(kTRUE);
        tlt.SetTextSize(0.55);
        std::ostringstream title;
        title << "Normalization grids: " << label << ", " << helicity
              << "   xB in [" << fmt_edge(xbk.xbmin_s) << ", " << fmt_edge(xbk.xbmax_s) << "]"
              << "   Ebeam=" << std::fixed << std::setprecision(1) << Ebeam
              << "   (log y)";
        tlt.DrawLatex(0.02, 0.20, title.str().c_str());

        c->cd();

        // Grid pad
        TPad *pgrid = new TPad((cname + "_grid").c_str(), "grid", 0.0, 0.0, 1.0, 0.92);
        pgrid->SetFillStyle(0);
        pgrid->SetFrameFillStyle(0);
        pgrid->Draw();
        pgrid->cd();

        pgrid->Divide(ncols, nrows, 0.001, 0.001);

        // Loop pads: row-major with q2 as rows (top to bottom), t as cols (left to right)
        for (int ir = 0; ir < nrows; ++ir) {
            for (int ic = 0; ic < ncols; ++ic) {
                const Q2BinKey &q2k = q2_bins[(size_t)ir];
                const TBinKey  &tk  = t_bins[(size_t)ic];

                KinKey kk;
                kk.xbmin_s = xbk.xbmin_s;
                kk.xbmax_s = xbk.xbmax_s;
                kk.q2min_s = q2k.q2min_s;
                kk.q2max_s = q2k.q2max_s;
                kk.tmin_s  = tk.tmin_s;
                kk.tmax_s  = tk.tmax_s;

                const int pad_idx = ir * ncols + ic + 1;
                TPad *pad = (TPad*)pgrid->cd(pad_idx);
                pad->SetLeftMargin(0.16);
                pad->SetRightMargin(0.07);
                pad->SetTopMargin(0.18);
                pad->SetBottomMargin(0.18);
                pad->SetGrid(1, 1);

                // NEW: log-y per subplot
                pad->SetLogy(1);

                // Subplot title line: Q2 and -t on one line
                const std::string q2tline =
                    "Q^{2} in [" + fmt_edge(q2k.q2min_s) + ", " + fmt_edge(q2k.q2max_s) + "], "
                    + "-t in [" + neg_edge_str(tk.tmax_s) + ", " + neg_edge_str(tk.tmin_s) + "]";

                std::map<KinKey, std::vector<PhiPoint> >::const_iterator itp = bin_points.find(kk);
                if (itp == bin_points.end()) {
                    draw_pad_text_only("No bin", q2tline, "");
                    continue;
                } //endif

                const std::vector<PhiPoint> &pts = itp->second;
                if (pts.empty()) {
                    draw_pad_text_only("No points", q2tline, "");
                    continue;
                } //endif

                // Keep graph alive via keep_graphs
                keep_graphs.emplace_back(std::make_unique<TGraphErrors>());
                TGraphErrors *gr = keep_graphs.back().get();

                gr->SetName(Form("gr_%s_r%d_c%d", cname.c_str(), ir, ic));
                gr->SetMarkerStyle(20);
                gr->SetMarkerSize(0.85);
                gr->SetMarkerColor(kBlack);
                gr->SetLineColor(kBlack);
                gr->SetLineWidth(1);

                int np_draw = 0;

                for (size_t i = 0; i < pts.size(); ++i) {
                    if (!std::isfinite(pts[i].phi_deg)) continue;
                    if (!std::isfinite(pts[i].xs) || pts[i].xs <= 0.0) continue;

                    const double ey = (std::isfinite(pts[i].xs_stat) && pts[i].xs_stat > 0.0 ? pts[i].xs_stat : 0.0);

                    gr->SetPoint(np_draw, pts[i].phi_deg, pts[i].xs);
                    gr->SetPointError(np_draw, 0.0, ey);
                    np_draw += 1;
                } //endfor

                if (np_draw <= 0) {
                    draw_pad_text_only("No usable points", q2tline, "");
                    continue;
                } //endif

                // Draw an explicit frame using the shared y-range
                TH1F *frame = (TH1F*)pad->DrawFrame(0.0, y_min_global, 360.0, y_max_global);
                frame->SetTitle("");
                frame->GetXaxis()->SetTitle("#phi (deg)");
                frame->GetYaxis()->SetTitle("xs");
                frame->GetXaxis()->SetNdivisions(505);
                frame->GetXaxis()->SetLabelSize(0.060);
                frame->GetYaxis()->SetLabelSize(0.060);
                frame->GetXaxis()->SetTitleSize(0.070);
                frame->GetYaxis()->SetTitleSize(0.070);
                frame->GetXaxis()->CenterTitle();
                frame->GetYaxis()->CenterTitle();

                // Subplot title (single line)
                {
                    TLatex tl;
                    tl.SetNDC(kTRUE);
                    tl.SetTextSize(0.055);
                    tl.DrawLatex(0.12, 0.92, q2tline.c_str());
                }

                // Draw points (with errors)
                gr->Draw("PE1 SAME");

                // Fit harmonics with iterative reduction
                int usedH = -1;
                double A = std::numeric_limits<double>::quiet_NaN();
                double dA = std::numeric_limits<double>::quiet_NaN();
                std::vector<double> pars;
                TMatrixDSym cov;
                std::string fit_msg;

                const bool ok_fit = fit_cos_series(pts, usedH, A, dA, pars, cov, fit_msg);

                // Bin centers (for BH/KM15 evaluation)
                double xb_c = std::numeric_limits<double>::quiet_NaN();
                double q2_c = std::numeric_limits<double>::quiet_NaN();
                double t_c  = std::numeric_limits<double>::quiet_NaN();

                std::map<KinKey, BinMeta>::const_iterator itm = bin_meta.find(kk);
                if (itm != bin_meta.end()) {
                    xb_c = itm->second.xb_c;
                    q2_c = itm->second.q2_c;
                    t_c  = itm->second.t_c;
                } //endif

                // Evaluate BH/KM15 at phi0
                double bh0 = std::numeric_limits<double>::quiet_NaN();
                double km0 = std::numeric_limits<double>::quiet_NaN();
                double bh_over_km0 = std::numeric_limits<double>::quiet_NaN();

                if (finite_pos(xb_c) && finite_pos(q2_c) && finite_pos(t_c)) {
                    bh0 = vgg_bh_only(xb_c, q2_c, t_c, phi0_deg, Ebeam);
                    km0 = km15_xs(xb_c, q2_c, t_c, phi0_deg, Ebeam, hel);
                    if (finite_pos(bh0) && finite_pos(km0)) {
                        bh_over_km0 = bh0 / km0;
                    } //endif
                } //endif

                // If fit ok, draw fitted function and annotate numbers
                if (ok_fit) {
                    // Keep fit function alive via keep_fits
                    const std::string fstr = build_cos_series_formula(usedH);
                    keep_fits.emplace_back(std::make_unique<TF1>(
                        Form("fdraw_%s_r%d_c%d", cname.c_str(), ir, ic),
                        fstr.c_str(),
                        0.0, 360.0
                    ));
                    TF1 *fdraw = keep_fits.back().get();

                    fdraw->SetNpx(720);
                    for (int ip = 0; ip <= usedH; ++ip) {
                        fdraw->SetParameter(ip, pars[(size_t)ip]);
                    } //endfor
                    fdraw->SetLineColor(kRed+1);
                    fdraw->SetLineWidth(2);
                    fdraw->Draw("SAME");

                    double xs_fit0 = std::numeric_limits<double>::quiet_NaN();
                    double dxs_fit0 = std::numeric_limits<double>::quiet_NaN();
                    const bool ok_eval = eval_fit_and_unc_at_phi0(phi0_deg, usedH, pars, cov, xs_fit0, dxs_fit0);

                    double sigma_over_bh0 = std::numeric_limits<double>::quiet_NaN();
                    if (finite_pos(bh0) && std::isfinite(dxs_fit0)) {
                        sigma_over_bh0 = dxs_fit0 / bh0;
                    } //endif

                    // Print row in summary table
                    if (finite_pos(xb_c) && finite_pos(q2_c) && finite_pos(t_c) &&
                        finite_pos(bh0) && ok_eval) {

                        std::cout
                            << std::setw(10) << xb_c
                            << std::setw(10) << q2_c
                            << std::setw(10) << (-t_c)
                            << std::setw(12) << A
                            << std::setw(12) << dA
                            << std::setw(12) << bh0
                            << std::setw(12) << (std::isfinite(bh_over_km0) ? bh_over_km0 : 0.0)
                            << std::setw(12) << xs_fit0
                            << std::setw(12) << dxs_fit0
                            << std::setw(14) << (std::isfinite(sigma_over_bh0) ? sigma_over_bh0 : 0.0)
                            << "\n";
                    } //endif

                    // Per-subplot legend: top-center, bordered, slightly smaller font
                    if (ok_eval && std::isfinite(bh_over_km0) && std::isfinite(dxs_fit0) && std::isfinite(sigma_over_bh0)) {
                        keep_legs.emplace_back(std::make_unique<TLegend>(0.24, 0.68, 0.76, 0.88));
                        TLegend *leg = keep_legs.back().get();

                        leg->SetBorderSize(1);
                        leg->SetFillStyle(1001);
                        leg->SetFillColor(kWhite);
                        leg->SetTextSize(0.045);

                        {
                            std::ostringstream s1;
                            s1 << "BH/KM15(#phi_{0}) = " << std::fixed << std::setprecision(3) << bh_over_km0;
                            leg->AddEntry((TObject*)0, s1.str().c_str(), "");
                        }
                        {
                            std::ostringstream s2;
                            s2 << "#sigma(#phi_{0}) = " << std::fixed << std::setprecision(5) << dxs_fit0;
                            leg->AddEntry((TObject*)0, s2.str().c_str(), "");
                        }
                        {
                            std::ostringstream s3;
                            s3 << "#sigma/BH(#phi_{0}) = " << std::fixed << std::setprecision(5) << sigma_over_bh0;
                            leg->AddEntry((TObject*)0, s3.str().c_str(), "");
                        }

                        leg->Draw();
                    } else {
                        TLatex tl;
                        tl.SetNDC(kTRUE);
                        tl.SetTextSize(0.050);
                        tl.DrawLatex(0.12, 0.80, "Eval failed");
                    } //endif
                } else {
                    // Fit failed: keep the data points visible (already drawn), annotate failure
                    TLatex tl;
                    tl.SetNDC(kTRUE);
                    tl.SetTextSize(0.055);
                    tl.DrawLatex(0.12, 0.80, "Fit failed");
                    if (!fit_msg.empty()) {
                        tl.SetTextSize(0.045);
                        tl.DrawLatex(0.12, 0.68, fit_msg.c_str());
                    } //endif
                } //endif
            } //endfor ic
        } //endfor ir

        // Save
        const std::string label_tag = sanitize_for_filename(label);
        const std::string hel_tag   = sanitize_for_filename(helicity);

        const std::string out_png =
            out_dir_grid
            + "/norm_grid_xb_"
            + fmt_edge_fname(xbk.xbmin_s) + "_"
            + fmt_edge_fname(xbk.xbmax_s) + "_"
            + label_tag + "_"
            + hel_tag + ".png";

        c->SaveAs(out_png.c_str());

        delete pgrid;
        delete ptitle;
        delete c;
    } //endfor xB bins

    std::cout << "------------------------------------------------------------\n";
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

        // ---------------------------------------------------------------------
        // Legacy selection: best phi (closest to 0/360) per (xB,Q2,t) bin
        // ---------------------------------------------------------------------
        std::map<std::string, BestPhiRow> best;

        // ---------------------------------------------------------------------
        // New: store all phi points per (xB,Q2,t) bin for harmonic fitting
        // ---------------------------------------------------------------------
        std::map<KinKey, std::vector<PhiPoint> > bin_points;
        std::map<KinKey, BinMeta> bin_meta;

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
            const double t_c    = std::atof(trim(unquote(fields[c_tavg])).c_str());
            const double phi    = std::atof(trim(unquote(fields[c_phiavg])).c_str());

            if (!finite_pos(xb_c) || !finite_pos(q2_c) || !finite_pos(t_c) || !std::isfinite(phi)) {
                continue;
            } //endif

            TripleCell xs = parse_tuple3(fields[c_xs]);
            if (!(xs.value > 0.0) || !std::isfinite(xs.value)) {
                continue;
            } //endif

            // -----------------------------------------------------------------
            // Legacy best-phi selection (unchanged)
            // -----------------------------------------------------------------
            {
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

            // -----------------------------------------------------------------
            // New: accumulate all points for this kinematic bin
            // -----------------------------------------------------------------
            {
                KinKey kk;
                kk.xbmin_s = xbmin_s;
                kk.xbmax_s = xbmax_s;
                kk.q2min_s = q2min_s;
                kk.q2max_s = q2max_s;
                kk.tmin_s  = tmin_s;
                kk.tmax_s  = tmax_s;

                PhiPoint pp;
                pp.phi_deg = phi;
                pp.xs = xs.value;
                pp.xs_stat = xs.stat;

                bin_points[kk].push_back(pp);

                // store centers once deterministically (same for all phi rows in that bin)
                if (bin_meta.find(kk) == bin_meta.end()) {
                    BinMeta bm;
                    bm.xb_c = xb_c;
                    bm.q2_c = q2_c;
                    bm.t_c  = t_c;
                    bin_meta[kk] = bm;
                } //endif
            }
        } //endfor rows

        if (best.empty()) {
            std::cerr << "[overall_norm] WARNING: no bins found with positive cross section values.\n";
            return true;
        } //endif

        std::cout << "[overall_norm] Unique kinematic bins found (legacy best-phi): " << best.size()
                  << " (from " << n_rows << " CSV data rows)\n\n";

        // ---------------------------------------------------------------------
        // Legacy printing and legacy plot-point creation (unchanged behavior)
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

        // Weighted mean accumulator for legacy xs/BH in BH/KM15 in [0.95,1.05]
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
                } //endif
            } //endif
            if (finite_pos(vgg)) {
                bh_over_vgg = bh / vgg;
            } //endif
            if (finite_pos(km)) {
                bh_over_km = bh / km;
            } //endif

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

            // Legacy plot-point selection: keep behavior identical:
            // require dist_to_edge > 0.0 (no explicit floor added)
            if (std::isfinite(br.dist_to_edge) && br.dist_to_edge > 0.0 &&
                std::isfinite(xs_over_bh) && xs_over_bh >= 0.0 &&
                std::isfinite(bh_over_km) && bh_over_km > 0.0) {

                PlotPoint p;
                p.d_edge_deg      = br.dist_to_edge;
                p.xs_over_bh      = xs_over_bh;
                p.xs_over_bh_err  = xs_over_bh_err;
                p.bh_over_km15    = bh_over_km;
                plot_pts.push_back(p);
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
        } //endfor best bins

        // ---------------------------------------------------------------------
        // Legacy plots (unchanged behavior/ranges)
        // ---------------------------------------------------------------------
        const std::string out_dir = "output/normalization_study";
        ensure_output_dir_or_throw(out_dir);
        make_normalization_plots_legacy(out_dir, label, helicity, plot_pts);

        // ---------------------------------------------------------------------
        // New grids (now: log-y + shared y-range per canvas)
        // ---------------------------------------------------------------------
        const std::string out_dir_grid = out_dir + "/grid";
        make_normalization_grids(out_dir_grid, label, helicity, Ebeam, hel, bin_points, bin_meta);

        // Legacy weighted mean summary (unchanged)
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