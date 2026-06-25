// pi0_contamination.cpp
// -----------------------------------------------------------------------------
// Pi0 contamination calculation using values already stored in the pass-2 CSV.
//
// This refactored implementation no longer loops over ROOT trees. Earlier
// modules have already filled:
//
//   1) normalized raw data yields for ep->epg and ep->eppi0
//   2) reconstructed current-corrected AAOGEN MC yields for ep->eppi0
//   3) reconstructed current-corrected AAOGEN background yields for
//      ep->eppi0->epg, i.e. generated pi0 events reconstructed as DVCS
//
// Therefore, for each period and each kinematic CSV row, the contamination is
// computed algebraically as:
//
//   c = N_mis * N_pi0_data / (N_pi0_rec_mc * N_dvcs_data)
//
// using topology-summed values. The CSV schema currently has one contamination
// ratio column per period, not one per topology, so the three standard topology
// columns are summed before computing the ratio.
// -----------------------------------------------------------------------------

#include "pi0_contamination.h"

#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TTree.h>
#include <TLeaf.h>
#include <TAxis.h>
#include <TBranch.h>
#include <TObjArray.h>

#include <algorithm>
#include <cmath>
#include <cctype>
#include <cstdio>
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
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

struct Triple {
    double v = 0.0;
    double stat = 0.0;
    double sys = 0.0;
};

struct CSV {
    std::vector<std::string> header;
    std::unordered_map<std::string, int> index;
    std::vector<std::vector<std::string>> rows;
};

struct RowBin {
    double xBmin = 0.0;
    double xBmax = 0.0;
    double Q2min = 0.0;
    double Q2max = 0.0;
    double tmin = 0.0;
    double tmax = 0.0;
    double phimin = 0.0;
    double phimax = 0.0;
    bool valid = false;
};

struct RowContam {
    double phi = 0.0;
    double value = 0.0;
    double stat = 0.0;
    double xBmin = 0.0;
    double xBmax = 0.0;
    double Q2min = 0.0;
    double Q2max = 0.0;
    double tmin = 0.0;
    double tmax = 0.0;
    bool valid = false;
};

static const std::vector<std::string> kPeriods = {
    "Fa18 Inb",
    "Fa18 Out",
    "Sp19 Inb",
    "Sp18 Inb",
    "Sp18 Out"
};

static const std::vector<std::string> kTopologies = {
    "(FD, FD)",
    "(CD, FD)",
    "(CD, FT)"
};

static void fatal(const std::string& msg) {
    throw std::runtime_error(msg);
}

static std::vector<std::string> split_csv_line(const std::string& line) {
    std::vector<std::string> out;
    std::string cur;
    cur.reserve(line.size());
    bool in_quotes = false;

    for (size_t i = 0; i < line.size(); ++i) {
        const char c = line[i];

        if (c == '"') {
            if (in_quotes && i + 1 < line.size() && line[i + 1] == '"') {
                cur.push_back('"');
                ++i;
            } else {
                in_quotes = !in_quotes;
            }
        } else if (c == ',' && !in_quotes) {
            out.push_back(cur);
            cur.clear();
        } else {
            cur.push_back(c);
        }
    }

    out.push_back(cur);
    return out;
}

static std::string join_csv_row(const std::vector<std::string>& fields) {
    std::ostringstream os;

    for (size_t i = 0; i < fields.size(); ++i) {
        const std::string& s = fields[i];
        const bool quote =
            s.find(',') != std::string::npos ||
            s.find('"') != std::string::npos ||
            s.find('\n') != std::string::npos ||
            s.find('\r') != std::string::npos;

        if (quote) {
            os << '"';
            for (char ch : s) {
                if (ch == '"') {
                    os << "\"\"";
                } else {
                    os << ch;
                }
            }
            os << '"';
        } else {
            os << s;
        }

        if (i + 1 < fields.size()) {
            os << ',';
        }
    }

    return os.str();
}

static void load_csv(const std::string& path, CSV& csv) {
    std::ifstream in(path);
    if (!in.is_open()) {
        fatal("[pi0_contamination] FATAL: cannot open CSV: " + path);
    }

    std::string line;
    if (!std::getline(in, line)) {
        fatal("[pi0_contamination] FATAL: empty CSV: " + path);
    }

    csv.header = split_csv_line(line);
    csv.index.clear();

    for (int i = 0; i < (int)csv.header.size(); ++i) {
        const std::string& name = csv.header[i];
        if (csv.index.find(name) != csv.index.end()) {
            fatal("[pi0_contamination] FATAL: duplicate CSV column: " + name);
        }
        csv.index[name] = i;
    }

    csv.rows.clear();
    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }

        std::vector<std::string> row = split_csv_line(line);
        if (row.size() < csv.header.size()) {
            row.resize(csv.header.size(), "");
        }
        if (row.size() != csv.header.size()) {
            fatal("[pi0_contamination] FATAL: CSV row width mismatch while reading: " + path);
        }
        csv.rows.push_back(std::move(row));
    }
}

static void write_csv_atomic(const std::string& path, const CSV& csv) {
    const std::string tmp = path + ".tmp";

    {
        std::ofstream out(tmp);
        if (!out.is_open()) {
            fatal("[pi0_contamination] FATAL: cannot write temp CSV: " + tmp);
        }

        out << join_csv_row(csv.header) << "\n";
        for (const auto& row : csv.rows) {
            out << join_csv_row(row) << "\n";
        }
    }

    (void)std::remove(path.c_str());
    if (std::rename(tmp.c_str(), path.c_str()) != 0) {
        fatal("[pi0_contamination] FATAL: atomic rename failed: " + tmp + " -> " + path);
    }
}

static int col_strict(const CSV& csv, const std::string& name) {
    auto it = csv.index.find(name);
    if (it == csv.index.end()) {
        fatal("[pi0_contamination] FATAL: missing CSV column: " + name);
    }
    return it->second;
}

static bool valid_cell(const std::string& s) {
    return s == "1" || s == "1.0" || s == "true" || s == "TRUE";
}

static double to_double_strict(const std::string& s, const std::string& context) {
    if (s.empty()) {
        fatal("[pi0_contamination] FATAL: empty numeric cell for " + context);
    }

    char* endp = nullptr;
    const double v = std::strtod(s.c_str(), &endp);

    if (endp == s.c_str() || *endp != '\0') {
        fatal("[pi0_contamination] FATAL: cannot parse numeric cell for " + context + ": '" + s + "'");
    }

    return v;
}

static bool parse_tuple_numbers(const std::string& s, std::vector<double>& vals) {
    vals.clear();

    std::string t;
    t.reserve(s.size());
    for (char c : s) {
        if (c == '(' || c == ')' || c == '"') {
            t.push_back(' ');
        } else {
            t.push_back(c);
        }
    }

    std::stringstream ss(t);
    std::string part;

    while (std::getline(ss, part, ',')) {
        char* endp = nullptr;
        const double v = std::strtod(part.c_str(), &endp);
        if (endp == part.c_str()) {
            return false;
        }
        vals.push_back(v);
    }

    return !vals.empty();
}

static Triple parse_yield_or_zero(const std::string& s,
                                  const std::string& context) {
    Triple out;
    out.v = 0.0;
    out.stat = 0.0;
    out.sys = 0.0;

    if (s.empty()) {
        return out;
    }

    std::vector<double> vals;
    if (!parse_tuple_numbers(s, vals) || vals.size() < 3) {
        fatal("[pi0_contamination] FATAL: expected yield triple (value,stat,sys) for " +
              context + ", got '" + s + "'.");
    }

    out.v = vals[0];
    out.stat = vals[1];
    out.sys = vals[2];

    if (!std::isfinite(out.v) || !std::isfinite(out.stat) || !std::isfinite(out.sys) ||
        out.stat < 0.0 || out.sys < 0.0) {
        fatal("[pi0_contamination] FATAL: invalid yield triple for " + context +
              ": '" + s + "'.");
    }

    return out;
}

static void add_triple_in_quadrature(Triple& total, const Triple& x) {
    total.v += x.v;
    total.stat = std::sqrt(total.stat * total.stat + x.stat * x.stat);
    total.sys = std::sqrt(total.sys * total.sys + x.sys * x.sys);
}

static std::string format_triple(const Triple& t) {
    std::ostringstream os;
    os << "("
       << std::fixed << std::setprecision(8) << t.v << ","
       << std::fixed << std::setprecision(8) << t.stat << ","
       << std::fixed << std::setprecision(8) << t.sys
       << ")";
    return os.str();
}

static std::vector<RowBin> load_row_bins(const CSV& csv) {
    const int c_xmin = col_strict(csv, "xBmin");
    const int c_xmax = col_strict(csv, "xBmax");
    const int c_qmin = col_strict(csv, "Q2min");
    const int c_qmax = col_strict(csv, "Q2max");
    const int c_tmin = col_strict(csv, "t_abs_min");
    const int c_tmax = col_strict(csv, "t_abs_max");
    const int c_pmin = col_strict(csv, "phimin");
    const int c_pmax = col_strict(csv, "phimax");
    const int c_valid = col_strict(csv, "valid bin");

    std::vector<RowBin> bins;
    bins.reserve(csv.rows.size());

    for (const auto& row : csv.rows) {
        RowBin b;
        b.xBmin = to_double_strict(row[c_xmin], "xBmin");
        b.xBmax = to_double_strict(row[c_xmax], "xBmax");
        b.Q2min = to_double_strict(row[c_qmin], "Q2min");
        b.Q2max = to_double_strict(row[c_qmax], "Q2max");
        b.tmin = to_double_strict(row[c_tmin], "t_abs_min");
        b.tmax = to_double_strict(row[c_tmax], "t_abs_max");
        b.phimin = to_double_strict(row[c_pmin], "phimin");
        b.phimax = to_double_strict(row[c_pmax], "phimax");
        b.valid = valid_cell(row[c_valid]);
        bins.push_back(b);
    }

    return bins;
}

static double phi_center(double pmin, double pmax) {
    if (pmax >= pmin) {
        return 0.5 * (pmin + pmax);
    }

    const double span = (360.0 - pmin) + pmax;
    double p = pmin + 0.5 * span;
    if (p >= 360.0) {
        p -= 360.0;
    }
    return p;
}

static std::string normalized_raw_yield_col(const std::string& channel,
                                            const std::string& topo,
                                            const std::string& period,
                                            const std::string& helicity) {
    return "normalized raw yield, " + channel + ", " + topo + ", exp, " + period + ", " + helicity;
}

static std::string rec_current_corrected_yield_col(const std::string& channel,
                                                   const std::string& topo,
                                                   const std::string& period) {
    return "reconstructed current corrected yield, " + channel + ", " + topo + ", mc, " + period;
}

static std::string contamination_col(const std::string& period) {
    return "contamination ratio, " + period;
}

static Triple sum_topology_data_yield(const CSV& csv,
                                      const std::vector<std::string>& row,
                                      const std::string& channel,
                                      const std::string& period) {
    Triple total;
    total.v = 0.0;
    total.stat = 0.0;
    total.sys = 0.0;

    for (const std::string& topo : kTopologies) {
        const std::string colname = normalized_raw_yield_col(channel, topo, period, "unpol");
        const int c = col_strict(csv, colname);
        add_triple_in_quadrature(total, parse_yield_or_zero(row[c], colname));
    }
    return total;
}

static Triple sum_topology_mc_yield(const CSV& csv,
                                    const std::vector<std::string>& row,
                                    const std::string& channel,
                                    const std::string& period) {
    Triple total;
    total.v = 0.0;
    total.stat = 0.0;
    total.sys = 0.0;

    for (const std::string& topo : kTopologies) {
        const std::string colname = rec_current_corrected_yield_col(channel, topo, period);
        const int c = col_strict(csv, colname);
        add_triple_in_quadrature(total, parse_yield_or_zero(row[c], colname));
    }
    return total;
}

static Triple compute_contamination(const Triple& Nmis,
                                    const Triple& Npi0_data,
                                    const Triple& Npi0_rec_mc,
                                    const Triple& Ndvcs_data) {
    Triple out;

    if (!(Nmis.v > 0.0) || !(Npi0_data.v > 0.0) ||
        !(Npi0_rec_mc.v > 0.0) || !(Ndvcs_data.v > 0.0)) {
        out.v = 0.0;
        out.stat = 0.0;
        out.sys = 0.0;
        return out;
    }

    out.v = (Nmis.v * Npi0_data.v) / (Npi0_rec_mc.v * Ndvcs_data.v);

    const double rel_var_stat =
        (Nmis.stat * Nmis.stat) / (Nmis.v * Nmis.v) +
        (Npi0_data.stat * Npi0_data.stat) / (Npi0_data.v * Npi0_data.v) +
        (Npi0_rec_mc.stat * Npi0_rec_mc.stat) / (Npi0_rec_mc.v * Npi0_rec_mc.v) +
        (Ndvcs_data.stat * Ndvcs_data.stat) / (Ndvcs_data.v * Ndvcs_data.v);

    out.stat = std::fabs(out.v) * std::sqrt(std::max(0.0, rel_var_stat));
    out.sys = 0.0;

    return out;
}

static std::string safe_dir_name(const std::string& s) {
    std::string out = s;
    for (char& c : out) {
        if (c == ' ') {
            c = '_';
        }
    }
    return out;
}

static std::string join_path(const std::string& a, const std::string& b) {
    if (a.empty()) {
        return b;
    }
    if (a.back() == '/') {
        return a + b;
    }
    return a + "/" + b;
}

static void mkdir_p(const std::string& path) {
    if (!path.empty()) {
        gSystem->mkdir(path.c_str(), true);
    }
}

static void plot_period_contamination(const std::string& period,
                                      const std::vector<RowBin>& bins,
                                      const std::vector<RowContam>& points,
                                      const std::string& output_root_dir) {
    const std::string out_dir = join_path(join_path(output_root_dir, "pi0_contamination_plots"), safe_dir_name(period));
    mkdir_p(out_dir);

    struct XBEdge {
        double a = 0.0;
        double b = 0.0;
    };

    std::vector<XBEdge> xbins;
    for (const RowBin& b : bins) {
        if (!b.valid) {
            continue;
        }
        XBEdge e{b.xBmin, b.xBmax};
        auto it = std::find_if(xbins.begin(), xbins.end(), [&](const XBEdge& z) {
            return z.a == e.a && z.b == e.b;
        });
        if (it == xbins.end()) {
            xbins.push_back(e);
        }
    }

    std::sort(xbins.begin(), xbins.end(), [](const XBEdge& x, const XBEdge& y) {
        if (x.a != y.a) {
            return x.a < y.a;
        }
        return x.b < y.b;
    });

    for (const XBEdge& xb : xbins) {
        std::vector<int> row_indices;
        std::vector<XBEdge> qbins;
        std::vector<XBEdge> tbins;

        for (int i = 0; i < (int)bins.size(); ++i) {
            const RowBin& b = bins[i];
            if (!b.valid) {
                continue;
            }
            if (!(b.xBmin == xb.a && b.xBmax == xb.b)) {
                continue;
            }

            row_indices.push_back(i);

            XBEdge q{b.Q2min, b.Q2max};
            XBEdge t{b.tmin, b.tmax};

            auto iq = std::find_if(qbins.begin(), qbins.end(), [&](const XBEdge& z) {
                return z.a == q.a && z.b == q.b;
            });
            if (iq == qbins.end()) {
                qbins.push_back(q);
            }

            auto it = std::find_if(tbins.begin(), tbins.end(), [&](const XBEdge& z) {
                return z.a == t.a && z.b == t.b;
            });
            if (it == tbins.end()) {
                tbins.push_back(t);
            }
        }

        if (row_indices.empty()) {
            continue;
        }

        auto sort_edge = [](std::vector<XBEdge>& v) {
            std::sort(v.begin(), v.end(), [](const XBEdge& x, const XBEdge& y) {
                if (x.a != y.a) {
                    return x.a < y.a;
                }
                return x.b < y.b;
            });
        };
        sort_edge(qbins);
        sort_edge(tbins);

        const int ncols = (int)qbins.size();
        const int nrows = (int)tbins.size();
        if (ncols <= 0 || nrows <= 0) {
            continue;
        }

        const int W = 300 * ncols + 160;
        const int H = 260 * nrows + 240;

        TCanvas c("c_pi0_contamination", "", W, H);
        TPad* top = new TPad("top", "top", 0.0, 0.90, 1.0, 1.0);
        TPad* grid = new TPad("grid", "grid", 0.0, 0.0, 1.0, 0.90);
        top->SetFillStyle(0);
        grid->SetFillStyle(0);
        top->Draw();
        grid->Draw();

        top->cd();
        TLatex title;
        title.SetNDC(true);
        title.SetTextFont(42);
        title.SetTextSize(0.45);
        std::ostringstream tss;
        tss << period << "   #pi_{0} contamination   x_{B} in ["
            << std::fixed << std::setprecision(2) << xb.a << ", " << xb.b << ")";
        title.DrawLatex(0.05, 0.35, tss.str().c_str());

        grid->cd();
        grid->Divide(ncols, nrows, 0.0, 0.0);

        auto find_edge = [](const std::vector<XBEdge>& v, double a, double b) {
            for (int i = 0; i < (int)v.size(); ++i) {
                if (v[i].a == a && v[i].b == b) {
                    return i;
                }
            }
            return -1;
        };

        for (int it = 0; it < nrows; ++it) {
            for (int iq = 0; iq < ncols; ++iq) {
                const int pad_idx = it * ncols + iq + 1;
                grid->cd(pad_idx);

                gPad->SetLeftMargin(0.160);
                gPad->SetRightMargin(0.07);
                gPad->SetTopMargin(0.22);
                gPad->SetBottomMargin(0.18);
                gPad->SetGrid(1, 1);
                gPad->SetTickx(1);
                gPad->SetTicky(1);

                TGraphErrors g;
                int ip = 0;
                double ymax = 0.0;

                for (int ridx : row_indices) {
                    const RowBin& rb = bins[ridx];
                    const int qidx = find_edge(qbins, rb.Q2min, rb.Q2max);
                    const int tidx = find_edge(tbins, rb.tmin, rb.tmax);
                    if (qidx != iq || tidx != it) {
                        continue;
                    }

                    const RowContam& p = points[ridx];
                    if (!p.valid) {
                        continue;
                    }

                    g.SetPoint(ip, p.phi, p.value);
                    g.SetPointError(ip, 0.0, p.stat);
                    ymax = std::max(ymax, p.value + p.stat);
                    ++ip;
                }

                if (!(ymax > 0.0)) {
                    ymax = 0.02;
                }

                TH1F* frame = (TH1F*)gPad->DrawFrame(0.0, 0.0, 360.0, 1.25 * ymax);
                frame->SetTitle("");
                frame->GetXaxis()->SetTitle("#phi (deg)");
                frame->GetYaxis()->SetTitle("#pi_{0} contamination");
                frame->GetXaxis()->CenterTitle(true);
                frame->GetYaxis()->CenterTitle(true);
                frame->GetXaxis()->SetNdivisions(505);
                frame->GetXaxis()->SetLabelSize(0.060);
                frame->GetYaxis()->SetLabelSize(0.060);
                frame->GetXaxis()->SetTitleSize(0.070);
                frame->GetYaxis()->SetTitleSize(0.070);
                frame->GetXaxis()->SetTitleOffset(1.05);
                frame->GetYaxis()->SetTitleOffset(1.10);

                g.SetMarkerStyle(20);
                g.SetMarkerColor(kBlack);
                g.SetLineColor(kBlack);
                g.SetLineWidth(1);
                g.Draw("PE SAME");

                TLatex txt;
                txt.SetNDC(true);
                txt.SetTextFont(42);
                txt.SetTextSize(0.060);
                std::ostringstream ss;
                ss << "Q^{2}=" << std::fixed << std::setprecision(2) << 0.5 * (qbins[iq].a + qbins[iq].b)
                   << "  |t|=" << std::fixed << std::setprecision(2) << 0.5 * (tbins[it].a + tbins[it].b);
                txt.DrawLatex(0.12, 0.83, ss.str().c_str());
            }
        }

        const int idx = (int)std::llround(xb.a * 1000.0);
        std::ostringstream fname;
        fname << out_dir << "/plot_pi0_contamination_" << safe_dir_name(period)
              << "_xB_" << idx << ".png";
        c.SaveAs(fname.str().c_str());

        delete top;
        delete grid;
    }
}


// -------------------- exclusivity-variable contamination diagnostics --------------------

struct DiagnosticVar {
    std::string name;
    std::string label;
    int nbins = 50;
    double xmin = 0.0;
    double xmax = 1.0;
};

struct DiagnosticBinRow {
    std::string period;
    std::string topology;
    std::string variable;
    int bin = 0;
    double xmin = 0.0;
    double xmax = 0.0;
    double xcenter = 0.0;
    double Ndvcs_data = 0.0;
    double Ndvcs_data_stat = 0.0;
    double Npi0_data = 0.0;
    double Npi0_data_stat = 0.0;
    double Npi0_rec_mc = 0.0;
    double Npi0_rec_mc_stat = 0.0;
    double Nmis_mc = 0.0;
    double Nmis_mc_stat = 0.0;
    double contamination = 0.0;       // fraction, not percent
    double contamination_stat = 0.0;  // fraction, not percent
    bool valid = false;
};

static const std::vector<DiagnosticVar> kDiagnosticVars = {
    {"Delta_phi",           "#Delta#phi (rad)",                                   100, 2.99, 3.29},
    {"theta_gamma_gamma",   "#theta_{#gamma#gamma} / #theta_{#pi^{0}#pi^{0}} (rad)", 100, 0.0, 2.0},
    {"pTmiss",              "p_{T}^{miss} (GeV)",                                 100, 0.0, 0.3},
    {"xF",                  "x_{F}",                                              100, -0.4, 0.2},
    {"Emiss2",              "E_{miss}^{2} (GeV^{2})",                             100, -1.0, 2.0},
    {"Mx2",                 "M_{x}^{2} (GeV^{2})",                                100, -0.03, 0.03},
    {"Mx2_1",               "M_{x1}^{2} (GeV^{2})",                               100, -1.5, 1.5},
    {"Mx2_2",               "M_{x2}^{2} (GeV^{2})",                               100, 0.0, 3.0}
};

static std::string to_lower_ascii(std::string s) {
    for (char& c : s) {
        c = (char)std::tolower((unsigned char)c);
    }
    return s;
}

static bool key_matches_period(const std::string& key, const std::string& period) {
    const std::string k = to_lower_ascii(key);

    if (period == "Fa18 Inb") return k.find("fa18") != std::string::npos && k.find("in")  != std::string::npos && k.find("out") == std::string::npos;
    if (period == "Fa18 Out") return k.find("fa18") != std::string::npos && k.find("out") != std::string::npos;
    if (period == "Sp18 Inb") return k.find("sp18") != std::string::npos && k.find("in")  != std::string::npos && k.find("out") == std::string::npos;
    if (period == "Sp18 Out") return k.find("sp18") != std::string::npos && k.find("out") != std::string::npos;
    if (period == "Sp19 Inb") return k.find("sp19") != std::string::npos && k.find("in")  != std::string::npos && k.find("out") == std::string::npos;

    return false;
}

static int topology_detector1(const std::string& topo) {
    if (topo == "(FD, FD)") return 1;
    if (topo == "(CD, FD)") return 2;
    if (topo == "(CD, FT)") return 2;
    return -999;
}

static int topology_detector2(const std::string& topo) {
    if (topo == "(FD, FD)") return 1;
    if (topo == "(CD, FD)") return 1;
    if (topo == "(CD, FT)") return 0;
    return -999;
}

static bool passes_topology_ids(int detector1, int detector2, const std::string& topo) {
    return detector1 == topology_detector1(topo) && detector2 == topology_detector2(topo);
}

static std::vector<std::string> diagnostic_variable_aliases(const std::string& requested_variable,
                                                            bool eppi0_exclusive_sample) {
    if (requested_variable == "theta_gamma_gamma") {
        if (eppi0_exclusive_sample) return {"theta_pi0_pi0", "theta_gamma_gamma"};
        return {"theta_gamma_gamma", "theta_pi0_pi0"};
    }

    if (requested_variable == "Delta_phi") {
        return {"Delta_phi"};
    }

    // Most trees use the nominal names below. The extra aliases are deliberately
    // harmless fallbacks for older skim variants. They also make the diagnostic
    // more robust if a channel-specific tree used a slightly different spelling.
    if (requested_variable == "Mx2") {
        return {"Mx2", "Mx2_epg", "Mx2_eg", "Mx2_epgamma", "Mx2_epi0", "Mx2_eppi0"};
    }
    if (requested_variable == "Mx2_1") {
        return {"Mx2_1", "Mx2_ep", "Mx2_x1", "Mx2_proton", "Mx2_p"};
    }
    if (requested_variable == "Mx2_2") {
        return {"Mx2_2", "Mx2_egamma", "Mx2_gamma", "Mx2_pi0", "Mx2_x2"};
    }

    return {requested_variable};
}

static std::string join_aliases_for_log(const std::vector<std::string>& aliases) {
    std::ostringstream os;
    for (size_t i = 0; i < aliases.size(); ++i) {
        if (i) os << "/";
        os << aliases[i];
    }
    return os.str();
}

static std::string safe_file_token(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    for (char c : s) {
        if (std::isalnum((unsigned char)c)) {
            out.push_back(c);
        } else {
            out.push_back('_');
        }
    }
    while (out.find("__") != std::string::npos) {
        size_t pos = out.find("__");
        out.replace(pos, 2, "_");
    }
    if (!out.empty() && out.front() == '_') out.erase(out.begin());
    if (!out.empty() && out.back() == '_') out.pop_back();
    return out.empty() ? "unnamed" : out;
}

struct TreeDiagLeaves {
    TTree* tree = nullptr;
    TLeaf* detector1 = nullptr;
    TLeaf* detector2 = nullptr;
    TLeaf* var = nullptr;
    TLeaf* Delta_phi = nullptr;
    TLeaf* p1_phi = nullptr;
    TLeaf* p2_phi = nullptr;
    std::string used_var_name;

    static std::string normalized_leaf_name(std::string s) {
        std::string out;
        out.reserve(s.size());
        for (char c : s) {
            if (c == '_' || c == '-' || c == ' ' || c == '.') continue;
            out.push_back((char)std::tolower((unsigned char)c));
        }
        return out;
    }

    static TLeaf* leaf_after_enabling(TTree* t, const std::string& name) {
        if (!t || name.empty()) return nullptr;

        // Previous stages sometimes disable branches to speed up tree scans.
        // Explicitly re-enable the branches needed by this diagnostic before
        // checking/filling them. Try both exact branch/leaf access and a final
        // case/underscore-insensitive leaf search for older tree variants.
        t->SetBranchStatus(name.c_str(), 1);

        if (TLeaf* leaf = t->GetLeaf(name.c_str())) return leaf;
        if (TBranch* branch = t->GetBranch(name.c_str())) {
            branch->SetStatus(1);
            if (TLeaf* leaf = branch->GetLeaf(name.c_str())) return leaf;
            TObjArray* leaves = branch->GetListOfLeaves();
            if (leaves && leaves->GetEntriesFast() == 1) {
                return dynamic_cast<TLeaf*>(leaves->At(0));
            }
        }

        const std::string want = normalized_leaf_name(name);
        TObjArray* leaves = t->GetListOfLeaves();
        if (!leaves) return nullptr;
        for (int i = 0; i < leaves->GetEntriesFast(); ++i) {
            TLeaf* leaf = dynamic_cast<TLeaf*>(leaves->At(i));
            if (!leaf) continue;
            const std::string lname = leaf->GetName() ? leaf->GetName() : "";
            if (normalized_leaf_name(lname) == want) {
                if (TBranch* branch = leaf->GetBranch()) {
                    branch->SetStatus(1);
                }
                return leaf;
            }
        }

        return nullptr;
    }

    bool init(TTree* t, const std::vector<std::string>& var_aliases) {
        tree = t;
        if (!tree) return false;

        detector1 = leaf_after_enabling(tree, "detector1");
        detector2 = leaf_after_enabling(tree, "detector2");
        Delta_phi = leaf_after_enabling(tree, "Delta_phi");
        p1_phi = leaf_after_enabling(tree, "p1_phi");
        p2_phi = leaf_after_enabling(tree, "p2_phi");

        for (const std::string& alias : var_aliases) {
            TLeaf* leaf = leaf_after_enabling(tree, alias);
            if (leaf) {
                var = leaf;
                used_var_name = alias;
                break;
            }
        }

        if (!detector1 || !detector2) return false;

        // Delta_phi can be stored directly or computed from p1_phi and p2_phi.
        if (var_aliases.size() == 1 && var_aliases.front() == "Delta_phi") {
            return (var != nullptr) || (Delta_phi != nullptr) || (p1_phi != nullptr && p2_phi != nullptr);
        }

        return var != nullptr;
    }

    static double wrapped_abs_delta_phi(double a, double b) {
        double d = a - b;
        constexpr double kPiLocal = 3.14159265358979323846;
        while (d >  kPiLocal) d -= 2.0 * kPiLocal;
        while (d < -kPiLocal) d += 2.0 * kPiLocal;
        return std::fabs(d);
    }

    bool get(double& value, int& d1, int& d2) const {
        if (!tree || !detector1 || !detector2) return false;
        d1 = (int)std::llround(detector1->GetValue());
        d2 = (int)std::llround(detector2->GetValue());

        if (var) {
            value = var->GetValue();
        } else if (Delta_phi) {
            value = Delta_phi->GetValue();
        } else if (p1_phi && p2_phi) {
            value = wrapped_abs_delta_phi(p1_phi->GetValue(), p2_phi->GetValue());
        } else {
            return false;
        }

        return std::isfinite(value);
    }
};

static void fill_diagnostic_hist_from_trees(const std::map<std::string, TTree*>& trees,
                                            const std::string& period,
                                            const std::string& topology,
                                            const DiagnosticVar& varcfg,
                                            bool eppi0_exclusive_sample,
                                            TH1D& hist,
                                            int& trees_used,
                                            long long& events_seen,
                                            long long& events_filled) {
    trees_used = 0;
    events_seen = 0;
    events_filled = 0;

    const std::vector<std::string> aliases = diagnostic_variable_aliases(varcfg.name, eppi0_exclusive_sample);
    static std::set<std::string> warned_missing;

    for (const auto& kv : trees) {
        if (!kv.second) continue;
        if (!key_matches_period(kv.first, period)) continue;

        TreeDiagLeaves leaves;
        if (!leaves.init(kv.second, aliases)) {
            const std::string warn_key = kv.first + "|" + varcfg.name + "|" + (eppi0_exclusive_sample ? "pi0" : "dvcs");
            if (warned_missing.insert(warn_key).second) {
                std::cout << "[pi0_contamination] diagnostic: skipping tree key='" << kv.first
                          << "' for variable " << varcfg.name
                          << " because required branch(es) are missing. Tried aliases="
                          << join_aliases_for_log(aliases) << "." << std::endl;
            }
            continue;
        }

        ++trees_used;
        const Long64_t n = kv.second->GetEntries();
        for (Long64_t i = 0; i < n; ++i) {
            kv.second->GetEntry(i);
            ++events_seen;

            double x = 0.0;
            int d1 = 0;
            int d2 = 0;
            if (!leaves.get(x, d1, d2)) continue;
            if (!passes_topology_ids(d1, d2, topology)) continue;
            if (x < varcfg.xmin || x >= varcfg.xmax) continue;

            hist.Fill(x);
            ++events_filled;
        }
    }
}

static Triple topology_data_yield_integrated(const CSV& csv,
                                             const std::vector<RowBin>& bins,
                                             const std::string& channel,
                                             const std::string& topo,
                                             const std::string& period) {
    Triple total;
    const std::string colname = normalized_raw_yield_col(channel, topo, period, "unpol");
    const int c = col_strict(csv, colname);

    for (int r = 0; r < (int)csv.rows.size(); ++r) {
        if (!bins[r].valid) continue;
        add_triple_in_quadrature(total, parse_yield_or_zero(csv.rows[r][c], colname));
    }
    return total;
}

static Triple topology_mc_yield_integrated(const CSV& csv,
                                           const std::vector<RowBin>& bins,
                                           const std::string& channel,
                                           const std::string& topo,
                                           const std::string& period) {
    Triple total;
    const std::string colname = rec_current_corrected_yield_col(channel, topo, period);
    const int c = col_strict(csv, colname);

    for (int r = 0; r < (int)csv.rows.size(); ++r) {
        if (!bins[r].valid) continue;
        add_triple_in_quadrature(total, parse_yield_or_zero(csv.rows[r][c], colname));
    }
    return total;
}

static void scale_hist_to_yield(TH1D& h, const Triple& target) {
    // Avoid TH1::Integral()/TH1::Scale() here. In this diagnostic pass we create
    // many short-lived stack histograms after long TTree scans. On some ROOT
    // builds TH1 can still have an internal fill buffer, and Integral() may call
    // BufferEmpty(), which has shown intermittent bus errors. The diagnostic only
    // needs ordinary bin sums, so do the scaling explicitly.
    double raw = 0.0;
    for (int ibin = 1; ibin <= h.GetNbinsX(); ++ibin) {
        raw += h.GetBinContent(ibin);
    }

    if (!(raw > 0.0) || !(target.v > 0.0)) return;

    const double scale = target.v / raw;
    for (int ibin = 1; ibin <= h.GetNbinsX(); ++ibin) {
        h.SetBinContent(ibin, h.GetBinContent(ibin) * scale);
        h.SetBinError(ibin, h.GetBinError(ibin) * scale);
    }
}

static DiagnosticBinRow compute_diagnostic_bin(const std::string& period,
                                               const std::string& topo,
                                               const std::string& var,
                                               int ibin,
                                               const TH1D& h_dvcs,
                                               const TH1D& h_pi0_data,
                                               const TH1D& h_pi0_rec,
                                               const TH1D& h_mis) {
    DiagnosticBinRow row;
    row.period = period;
    row.topology = topo;
    row.variable = var;
    row.bin = ibin;
    row.xmin = h_dvcs.GetXaxis()->GetBinLowEdge(ibin);
    row.xmax = h_dvcs.GetXaxis()->GetBinUpEdge(ibin);
    row.xcenter = h_dvcs.GetXaxis()->GetBinCenter(ibin);

    row.Ndvcs_data = h_dvcs.GetBinContent(ibin);
    row.Ndvcs_data_stat = h_dvcs.GetBinError(ibin);
    row.Npi0_data = h_pi0_data.GetBinContent(ibin);
    row.Npi0_data_stat = h_pi0_data.GetBinError(ibin);
    row.Npi0_rec_mc = h_pi0_rec.GetBinContent(ibin);
    row.Npi0_rec_mc_stat = h_pi0_rec.GetBinError(ibin);
    row.Nmis_mc = h_mis.GetBinContent(ibin);
    row.Nmis_mc_stat = h_mis.GetBinError(ibin);

    if (!(row.Ndvcs_data > 0.0) || !(row.Npi0_data > 0.0) ||
        !(row.Npi0_rec_mc > 0.0) || !(row.Nmis_mc > 0.0)) {
        row.valid = false;
        return row;
    }

    row.contamination = (row.Nmis_mc * row.Npi0_data) /
                        (row.Npi0_rec_mc * row.Ndvcs_data);

    const double rel_var =
        (row.Nmis_mc_stat * row.Nmis_mc_stat) / (row.Nmis_mc * row.Nmis_mc) +
        (row.Npi0_data_stat * row.Npi0_data_stat) / (row.Npi0_data * row.Npi0_data) +
        (row.Npi0_rec_mc_stat * row.Npi0_rec_mc_stat) / (row.Npi0_rec_mc * row.Npi0_rec_mc) +
        (row.Ndvcs_data_stat * row.Ndvcs_data_stat) / (row.Ndvcs_data * row.Ndvcs_data);

    row.contamination_stat = std::fabs(row.contamination) * std::sqrt(std::max(0.0, rel_var));
    row.valid = std::isfinite(row.contamination) && std::isfinite(row.contamination_stat);
    return row;
}

static void write_diagnostic_csv(const std::string& path,
                                 const std::vector<DiagnosticBinRow>& rows) {
    std::ofstream out(path);
    if (!out.is_open()) {
        fatal("[pi0_contamination] FATAL: cannot write diagnostic CSV: " + path);
    }

    out << "period,topology,variable,bin,xmin,xmax,xcenter,"
        << "Ndvcs_data,Ndvcs_data_stat,Npi0_data,Npi0_data_stat,"
        << "Npi0_rec_mc,Npi0_rec_mc_stat,Nmis_mc,Nmis_mc_stat,"
        << "contamination_fraction,contamination_fraction_stat,"
        << "contamination_percent,contamination_percent_stat,valid\n";

    out << std::setprecision(12);
    for (const DiagnosticBinRow& r : rows) {
        out << r.period << ',' << r.topology << ',' << r.variable << ',' << r.bin << ','
            << r.xmin << ',' << r.xmax << ',' << r.xcenter << ','
            << r.Ndvcs_data << ',' << r.Ndvcs_data_stat << ','
            << r.Npi0_data << ',' << r.Npi0_data_stat << ','
            << r.Npi0_rec_mc << ',' << r.Npi0_rec_mc_stat << ','
            << r.Nmis_mc << ',' << r.Nmis_mc_stat << ','
            << r.contamination << ',' << r.contamination_stat << ','
            << 100.0 * r.contamination << ',' << 100.0 * r.contamination_stat << ','
            << (r.valid ? 1 : 0) << '\n';
    }
}

static void draw_topology_graph(const std::vector<DiagnosticBinRow>& rows,
                                int color,
                                int marker,
                                const char* draw_opt,
                                TLegend& leg,
                                const std::string& topo) {
    TGraphErrors* g = new TGraphErrors();
    g->SetName(("g_" + safe_file_token(topo)).c_str());

    int ip = 0;
    for (const DiagnosticBinRow& r : rows) {
        if (!r.valid) continue;
        g->SetPoint(ip, r.xcenter, 100.0 * r.contamination);
        g->SetPointError(ip, 0.5 * (r.xmax - r.xmin), 100.0 * r.contamination_stat);
        ++ip;
    }

    g->SetMarkerStyle(marker);
    g->SetMarkerSize(0.75);
    g->SetMarkerColor(color);
    g->SetLineColor(color);
    g->SetLineWidth(1);

    if (ip > 0) {
        g->Draw(draw_opt);
        leg.AddEntry(g, topo.c_str(), "pe");
    }
}

static void plot_diagnostic_overlay(const std::string& out_png,
                                    const std::string& period,
                                    const DiagnosticVar& varcfg,
                                    const std::map<std::string, std::vector<DiagnosticBinRow>>& rows_by_topology) {
    TCanvas c("c_pi0_contamination_vs_excl_overlay", "", 1100, 800);
    c.SetLeftMargin(0.12);
    c.SetRightMargin(0.05);
    c.SetTopMargin(0.10);
    c.SetBottomMargin(0.13);
    c.SetGrid(1, 1);
    c.SetTickx(1);
    c.SetTicky(1);

    TH1F* frame = (TH1F*)gPad->DrawFrame(varcfg.xmin, 0.0, varcfg.xmax, 100.0);
    frame->SetTitle("");
    frame->GetXaxis()->SetTitle(varcfg.label.c_str());
    frame->GetYaxis()->SetTitle("predicted #pi_{0} contamination (%)");
    frame->GetXaxis()->CenterTitle(true);
    frame->GetYaxis()->CenterTitle(true);
    frame->GetXaxis()->SetTitleSize(0.045);
    frame->GetYaxis()->SetTitleSize(0.045);
    frame->GetXaxis()->SetLabelSize(0.040);
    frame->GetYaxis()->SetLabelSize(0.040);
    frame->GetYaxis()->SetTitleOffset(1.25);

    TLegend leg(0.68, 0.70, 0.93, 0.88);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.SetTextFont(42);
    leg.SetTextSize(0.035);

    const std::vector<int> colors = {kBlack, kBlue + 1, kRed + 1};
    const std::vector<int> markers = {20, 21, 22};

    int itopo = 0;
    for (const std::string& topo : kTopologies) {
        auto it = rows_by_topology.find(topo);
        if (it == rows_by_topology.end()) continue;
        draw_topology_graph(it->second,
                            colors[itopo % (int)colors.size()],
                            markers[itopo % (int)markers.size()],
                            "PE SAME",
                            leg,
                            topo);
        ++itopo;
    }

    leg.Draw();

    TLatex title;
    title.SetNDC(true);
    title.SetTextFont(42);
    title.SetTextSize(0.040);
    std::ostringstream ss;
    ss << period << "   predicted #pi_{0} contamination vs " << varcfg.name;
    title.DrawLatex(0.12, 0.935, ss.str().c_str());

    c.SaveAs(out_png.c_str());
}

static void make_exclusivity_variable_diagnostics(
    const CSV& csv,
    const std::vector<RowBin>& bins,
    const std::map<std::string, TTree*>& dvcsDataTrees,
    const std::map<std::string, TTree*>& eppi0DataTrees,
    const std::map<std::string, TTree*>& eppi0RecMcTrees,
    const std::map<std::string, TTree*>& eppi0BkgTrees,
    const std::string& output_root_dir) {

    (void)csv;
    (void)bins;

    const std::string diag_dir = join_path(join_path(output_root_dir, "pi0_contamination_plots"),
                                           "exclusivity_variable_diagnostics");
    mkdir_p(diag_dir);

    std::vector<DiagnosticBinRow> all_rows;

    for (const std::string& period : kPeriods) {
        for (const DiagnosticVar& varcfg : kDiagnosticVars) {
            std::map<std::string, std::vector<DiagnosticBinRow>> rows_by_topology;

            for (const std::string& topo : kTopologies) {
                const std::string hbase = safe_file_token(period + "_" + topo + "_" + varcfg.name);

                TH1D h_dvcs(("h_dvcs_" + hbase).c_str(), "", varcfg.nbins, varcfg.xmin, varcfg.xmax);
                TH1D h_pi0_data(("h_pi0_data_" + hbase).c_str(), "", varcfg.nbins, varcfg.xmin, varcfg.xmax);
                TH1D h_pi0_rec(("h_pi0_rec_" + hbase).c_str(), "", varcfg.nbins, varcfg.xmin, varcfg.xmax);
                TH1D h_mis(("h_mis_" + hbase).c_str(), "", varcfg.nbins, varcfg.xmin, varcfg.xmax);
                h_dvcs.SetDirectory(nullptr);
                h_pi0_data.SetDirectory(nullptr);
                h_pi0_rec.SetDirectory(nullptr);
                h_mis.SetDirectory(nullptr);
                h_dvcs.SetBuffer(0);
                h_pi0_data.SetBuffer(0);
                h_pi0_rec.SetBuffer(0);
                h_mis.SetBuffer(0);
                h_dvcs.Sumw2();
                h_pi0_data.Sumw2();
                h_pi0_rec.Sumw2();
                h_mis.Sumw2();

                int nt_dvcs = 0, nt_pi0_data = 0, nt_pi0_rec = 0, nt_mis = 0;
                long long ns_dvcs = 0, ns_pi0_data = 0, ns_pi0_rec = 0, ns_mis = 0;
                long long nf_dvcs = 0, nf_pi0_data = 0, nf_pi0_rec = 0, nf_mis = 0;

                fill_diagnostic_hist_from_trees(dvcsDataTrees, period, topo, varcfg,
                                                false, h_dvcs, nt_dvcs, ns_dvcs, nf_dvcs);
                fill_diagnostic_hist_from_trees(eppi0DataTrees, period, topo, varcfg,
                                                true, h_pi0_data, nt_pi0_data, ns_pi0_data, nf_pi0_data);
                fill_diagnostic_hist_from_trees(eppi0RecMcTrees, period, topo, varcfg,
                                                true, h_pi0_rec, nt_pi0_rec, ns_pi0_rec, nf_pi0_rec);
                fill_diagnostic_hist_from_trees(eppi0BkgTrees, period, topo, varcfg,
                                                false, h_mis, nt_mis, ns_mis, nf_mis);

                const bool have_all_trees = (nt_dvcs > 0 && nt_pi0_data > 0 && nt_pi0_rec > 0 && nt_mis > 0);

                if (!have_all_trees) {
                    std::cout << "[pi0_contamination] diagnostic: insufficient trees for "
                              << period << " " << topo << " " << varcfg.name
                              << " trees(dvcs,pi0data,pi0rec,mis)=("
                              << nt_dvcs << "," << nt_pi0_data << "," << nt_pi0_rec << "," << nt_mis << ")."
                              << std::endl;
                }

                // Do not renormalize these diagnostic histograms to CSV-integrated
                // totals. This study is the same pi0-contamination formula evaluated
                // directly as a function of the plotted exclusivity variable:
                //
                //   c(x) = Nmis_mc(x) * Npi0_data(x)
                //        / [Npi0_rec_mc(x) * Ndvcs_data(x)]
                //
                // i.e. the same structure used in the main calculation, but binned
                // in Delta_phi, pTmiss, Mx2, etc. rather than in xB, Q2, |t|, phi.

                std::vector<DiagnosticBinRow> plot_rows;
                plot_rows.reserve(varcfg.nbins);
                for (int ibin = 1; ibin <= varcfg.nbins; ++ibin) {
                    DiagnosticBinRow row = compute_diagnostic_bin(period, topo, varcfg.name, ibin,
                                                                  h_dvcs, h_pi0_data, h_pi0_rec, h_mis);
                    if (!have_all_trees) {
                        row.valid = false;
                    }
                    all_rows.push_back(row);
                    plot_rows.push_back(row);
                }

                rows_by_topology[topo] = plot_rows;

                std::cout << "[pi0_contamination] diagnostic " << period << " " << topo << " " << varcfg.name
                          << ": seen(dvcs,pi0data,pi0rec,mis)=("
                          << ns_dvcs << "," << ns_pi0_data << "," << ns_pi0_rec << "," << ns_mis
                          << ") filled(dvcs,pi0data,pi0rec,mis)=("
                          << nf_dvcs << "," << nf_pi0_data << "," << nf_pi0_rec << "," << nf_mis
                          << ")" << std::endl;
            }

            const std::string out_png = join_path(diag_dir,
                "pi0_contamination_vs_" + safe_file_token(varcfg.name) + "_" +
                safe_file_token(period) + "_topologies.png");
            plot_diagnostic_overlay(out_png, period, varcfg, rows_by_topology);
        }
    }

    const std::string csv_path = join_path(diag_dir, "pi0_contamination_vs_exclusivity_variables.csv");
    write_diagnostic_csv(csv_path, all_rows);

    std::cout << "[pi0_contamination] Wrote exclusivity-variable contamination diagnostics under: "
              << diag_dir << std::endl;
}

static void compute_period(CSV& csv,
                           const std::vector<RowBin>& bins,
                           const std::string& period,
                           const std::string& output_root_dir) {
    const int c_contam = col_strict(csv, contamination_col(period));

    std::vector<RowContam> points;
    points.resize(csv.rows.size());

    Triple sum_Ndvcs;
    Triple sum_Nmis;
    Triple sum_Npi0_data;
    Triple sum_Npi0_rec;
    sum_Ndvcs.v = sum_Ndvcs.stat = sum_Ndvcs.sys = 0.0;
    sum_Nmis.v = sum_Nmis.stat = sum_Nmis.sys = 0.0;
    sum_Npi0_data.v = sum_Npi0_data.stat = sum_Npi0_data.sys = 0.0;
    sum_Npi0_rec.v = sum_Npi0_rec.stat = sum_Npi0_rec.sys = 0.0;

    for (int r = 0; r < (int)csv.rows.size(); ++r) {
        const RowBin& b = bins[r];

        if (!b.valid) {
            csv.rows[r][c_contam].clear();
            continue;
        }

        const std::vector<std::string>& row = csv.rows[r];

        const Triple Ndvcs = sum_topology_data_yield(csv, row, "ep->epg", period);
        const Triple Npi0_data = sum_topology_data_yield(csv, row, "ep->eppi0", period);
        const Triple Npi0_rec = sum_topology_mc_yield(csv, row, "ep->eppi0", period);
        const Triple Nmis = sum_topology_mc_yield(csv, row, "ep->eppi0->epg", period);

        const Triple tr = compute_contamination(Nmis, Npi0_data, Npi0_rec, Ndvcs);
        csv.rows[r][c_contam] = format_triple(tr);

        RowContam p;
        p.phi = phi_center(b.phimin, b.phimax);
        p.value = tr.v;
        p.stat = tr.stat;
        p.xBmin = b.xBmin;
        p.xBmax = b.xBmax;
        p.Q2min = b.Q2min;
        p.Q2max = b.Q2max;
        p.tmin = b.tmin;
        p.tmax = b.tmax;
        p.valid = true;
        points[r] = p;

        add_triple_in_quadrature(sum_Ndvcs, Ndvcs);
        add_triple_in_quadrature(sum_Npi0_data, Npi0_data);
        add_triple_in_quadrature(sum_Npi0_rec, Npi0_rec);
        add_triple_in_quadrature(sum_Nmis, Nmis);
    }

    const Triple integrated = compute_contamination(sum_Nmis, sum_Npi0_data, sum_Npi0_rec, sum_Ndvcs);

    std::cout << "[pi0_contamination] " << period
              << " totals: Ndvcs_norm=" << std::setprecision(10) << sum_Ndvcs.v
              << " Nmis_corr_mc=" << sum_Nmis.v
              << " Npi0_norm=" << sum_Npi0_data.v
              << " Npi0_rec_corr_mc=" << sum_Npi0_rec.v
              << " integrated_contam=" << integrated.v
              << " +/- " << integrated.stat
              << std::endl;

    plot_period_contamination(period, bins, points, output_root_dir);
}

} // namespace

bool compute_pi0_contamination_overall(
    const std::map<std::string, TTree*> &dvcsDataTrees,
    const std::map<std::string, TTree*> &eppi0DataTrees,
    const std::map<std::string, TTree*> &eppi0RecMcTrees,
    const std::map<std::string, TTree*> &eppi0BkgTrees,
    const std::string &combined_cuts_json,
    const std::string &csv_main,
    const std::string &output_root_dir,
    int max_workers) {
    (void)combined_cuts_json;
    (void)max_workers;

    try {
        TH1::AddDirectory(kFALSE);
        gStyle->SetOptStat(0);

        CSV csv;
        load_csv(csv_main, csv);
        const std::vector<RowBin> bins = load_row_bins(csv);

        if (bins.size() != csv.rows.size()) {
            fatal("[pi0_contamination] FATAL: internal row/bin size mismatch.");
        }

        for (const std::string& period : kPeriods) {
            compute_period(csv, bins, period, output_root_dir);
        }

        make_exclusivity_variable_diagnostics(csv,
                                              bins,
                                              dvcsDataTrees,
                                              eppi0DataTrees,
                                              eppi0RecMcTrees,
                                              eppi0BkgTrees,
                                              output_root_dir);

        write_csv_atomic(csv_main, csv);

        std::cout << "[pi0_contamination] Updated contamination columns in: "
                  << csv_main << std::endl;
        std::cout << "[pi0_contamination] Plots written under: "
                  << join_path(output_root_dir, "pi0_contamination_plots") << std::endl;

        return true;
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return false;
    }
}
