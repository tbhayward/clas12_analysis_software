#include "overall_normalization_study.h"

#include "model_predictions.h"

// ROOT includes
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TProfile2D.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TH1F.h>
#include <TLine.h>

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
    double t_c;
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

static std::string sanitize_filename_component(const std::string &s) {
    std::string out;
    out.reserve(s.size());
    for (size_t i = 0; i < s.size(); ++i) {
        const unsigned char c = (unsigned char)s[i];
        if (std::isalnum(c)) {
            out.push_back((char)c);
        } else if (c == ' ' || c == '-' || c == '_' ) {
            out.push_back('_');
        } else {
            // drop other punctuation
        }
    }
    if (out.empty()) out = "NA";
    return out;
}

static void ensure_dir_or_throw(const std::string &dir) {
    int rc = gSystem->mkdir(dir.c_str(), kTRUE);
    if (rc != 0) {
        throw std::runtime_error("Failed to create output directory: " + dir);
    }
}

struct NormPoint {
    double d_edge;
    double bh_over_km15;
    double xs_over_bh;
    double xs_over_bh_err;
};

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

        std::vector<NormPoint> pts;
        pts.reserve(best.size());

        double min_d = 1e300, max_d = -1e300;
        double min_bhk = 1e300, max_bhk = -1e300;
        double min_xsb = 1e300, max_xsb = -1e300;

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
            double bh_over_vgg = 0.0;
            double bh_over_km  = 0.0;
            double xs_over_bh_err = 0.0;

            if (finite_pos(bh)) {
                xs_over_bh = br.xs / bh;
                if (br.xs_stat > 0.0 && std::isfinite(br.xs_stat)) {
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

            if (std::isfinite(br.dist_to_edge) && std::isfinite(bh_over_km) && std::isfinite(xs_over_bh)) {
                NormPoint p;
                p.d_edge = br.dist_to_edge;
                p.bh_over_km15 = bh_over_km;
                p.xs_over_bh = xs_over_bh;
                p.xs_over_bh_err = xs_over_bh_err;

                pts.push_back(p);

                if (p.d_edge < min_d) min_d = p.d_edge;
                if (p.d_edge > max_d) max_d = p.d_edge;
                if (p.bh_over_km15 < min_bhk) min_bhk = p.bh_over_km15;
                if (p.bh_over_km15 > max_bhk) max_bhk = p.bh_over_km15;
                if (p.xs_over_bh < min_xsb) min_xsb = p.xs_over_bh;
                if (p.xs_over_bh > max_xsb) max_xsb = p.xs_over_bh;
            }
        }

        // ---------------- Plotting block ----------------
        {
            const std::string outdir = "output/normalization_study";
            ensure_dir_or_throw(outdir);

            const std::string safe_label = sanitize_filename_component(label);
            const std::string safe_hel   = sanitize_filename_component(helicity);

            if (!(std::isfinite(min_d) && std::isfinite(max_d) && max_d > min_d)) {
                min_d = 0.0;
                max_d = 180.0;
            }
            if (!(std::isfinite(min_bhk) && std::isfinite(max_bhk) && max_bhk > min_bhk)) {
                min_bhk = 0.0;
                max_bhk = 1.5;
            }
            if (!(std::isfinite(min_xsb) && std::isfinite(max_xsb) && max_xsb > min_xsb)) {
                min_xsb = 0.0;
                max_xsb = 3.0;
            }

            const double d_pad   = 0.05 * (max_d - min_d);
            const double bhk_pad = 0.10 * (max_bhk - min_bhk);

            const double d_lo   = std::max(0.0, min_d - d_pad);
            const double d_hi   = max_d + d_pad;
            const double bhk_lo = std::max(0.0, min_bhk - bhk_pad);
            const double bhk_hi = max_bhk + bhk_pad;

            gStyle->SetOptStat(0);
            gStyle->SetTitleFontSize(0.045);

            // -------- Plot A: 2D heatmap (BH/KM15 vs d_edge), color = xs/BH --------
            {
                const int nx = 60;
                const int ny = 60;

                TCanvas *c = new TCanvas("c_norm_heatmap", "c_norm_heatmap", 1000, 800);
                c->SetLeftMargin(0.12);
                c->SetRightMargin(0.14);
                c->SetBottomMargin(0.12);
                c->SetTopMargin(0.10);

                std::ostringstream title;
                title << "Normalization study: xs/BH color; label=" << label << "; hel=" << helicity;

                TProfile2D *p = new TProfile2D(
                    "p_norm_xsOverBH",
                    title.str().c_str(),
                    nx, d_lo, d_hi,
                    ny, bhk_lo, bhk_hi
                );
                p->GetXaxis()->SetTitle("d_edge (deg)");
                p->GetYaxis()->SetTitle("BH/KM15");
                p->GetZaxis()->SetTitle("xs/BH");

                for (size_t i = 0; i < pts.size(); ++i) {
                    p->Fill(pts[i].d_edge, pts[i].bh_over_km15, pts[i].xs_over_bh);
                }

                p->Draw("COLZ");

                TLatex latex;
                latex.SetNDC();
                latex.SetTextSize(0.035);
                latex.DrawLatex(0.12, 0.93, (std::string("label: ") + label + "    hel: " + helicity).c_str());
                // N=... removed

                const std::string out = outdir + "/norm_heatmap_xsOverBH__label_" + safe_label + "__hel_" + safe_hel + ".png";
                c->SaveAs(out.c_str());

                delete p;
                delete c;
            }

            // -------- Plot B: 1D scatter xs/BH vs d_edge, grouped by BH/KM15 ranges --------
            {
                const double edges[] = {0.0, 0.70, 0.85, 1.00, 1.15, 1e9};
                const int ncat = 5;

                TGraphErrors *gr[ncat];
                for (int i = 0; i < ncat; ++i) {
                    std::ostringstream name;
                    name << "gr_norm_cat_" << i;
                    gr[i] = new TGraphErrors();
                    gr[i]->SetName(name.str().c_str());
                }

                auto cat_index = [&](double bhk) -> int {
                    for (int i = 0; i < ncat; ++i) {
                        if (bhk >= edges[i] && bhk < edges[i + 1]) return i;
                    }
                    return ncat - 1;
                };

                for (size_t i = 0; i < pts.size(); ++i) {
                    const int ci = cat_index(pts[i].bh_over_km15);
                    const int n  = gr[ci]->GetN();
                    gr[ci]->SetPoint(n, pts[i].d_edge, pts[i].xs_over_bh);
                    gr[ci]->SetPointError(n, 0.0, pts[i].xs_over_bh_err);
                }

                // Force y-range [0, 3] as requested
                const double y_lo = 0.0;
                const double y_hi = 3.0;

                TCanvas *c = new TCanvas("c_norm_scatter", "c_norm_scatter", 1100, 750);
                c->SetLeftMargin(0.12);
                c->SetRightMargin(0.05);
                c->SetBottomMargin(0.12);
                c->SetTopMargin(0.10);

                std::ostringstream frame_title;
                frame_title << "Normalization study: xs/BH vs d_edge; label=" << label << "; hel=" << helicity;

                TH1F *frame = new TH1F("frame_norm_scatter", frame_title.str().c_str(), 100, d_lo, d_hi);
                frame->GetXaxis()->SetTitle("d_edge (deg)");
                frame->GetYaxis()->SetTitle("xs/BH");
                frame->SetMinimum(y_lo);
                frame->SetMaximum(y_hi);
                frame->Draw("AXIS");

                // Add dashed gray line at y = 1
                TLine *line = new TLine(d_lo, 1.0, d_hi, 1.0);
                line->SetLineStyle(2);
                line->SetLineWidth(2);
                line->SetLineColor(16); // gray-ish (ROOT palette index)
                line->Draw("SAME");

                const int mstyles[ncat] = {20, 21, 22, 23, 24};
                const int mcolors[ncat] = {4, 38, 1, 2, 6};

                TLegend *leg = new TLegend(0.62, 0.70, 0.92, 0.90);
                leg->SetFillStyle(1001);
                leg->SetFillColor(0);
                leg->SetBorderSize(1);
                leg->SetTextSize(0.030);

                for (int i = 0; i < ncat; ++i) {
                    gr[i]->SetMarkerStyle(mstyles[i]);
                    gr[i]->SetMarkerSize(1.0);
                    gr[i]->SetMarkerColor(mcolors[i]);
                    gr[i]->SetLineColor(mcolors[i]);

                    if (gr[i]->GetN() > 0) {
                        gr[i]->Draw("P SAME");
                    }

                    std::ostringstream lab;
                    if (i < ncat - 1 && edges[i + 1] < 1e8) {
                        lab << "BH/KM15 in [" << std::fixed << std::setprecision(2) << edges[i]
                            << ", " << edges[i + 1] << ")";
                    } else {
                        lab << "BH/KM15 >= " << std::fixed << std::setprecision(2) << edges[i];
                    }
                    leg->AddEntry(gr[i], lab.str().c_str(), "p");
                }

                leg->Draw();

                TLatex latex;
                latex.SetNDC();
                latex.SetTextSize(0.035);
                latex.DrawLatex(0.12, 0.93, (std::string("label: ") + label + "    hel: " + helicity).c_str());
                // N=... removed

                const std::string out = outdir + "/norm_scatter_xsOverBH_vs_dedge__label_" + safe_label + "__hel_" + safe_hel + ".png";
                c->SaveAs(out.c_str());

                delete leg;
                delete line;
                delete frame;
                for (int i = 0; i < ncat; ++i) {
                    delete gr[i];
                }
                delete c;
            }
        }
        // ---------------- End plotting block ----------------

        std::cout << "\n";
        std::cout << "[overall_norm] Done.\n";
        std::cout << "============================================================\n";
        return true;

    } catch (const std::exception &e) {
        std::cerr << "[overall_norm] FATAL: " << e.what() << "\n";
        return false;
    }
}