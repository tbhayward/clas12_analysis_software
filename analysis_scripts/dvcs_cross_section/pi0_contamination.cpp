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
#include <TAxis.h>

#include <algorithm>
#include <cmath>
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

static double to_double_or_zero(const std::string& s) {
    if (s.empty()) {
        return 0.0;
    }

    char* endp = nullptr;
    const double v = std::strtod(s.c_str(), &endp);

    if (endp == s.c_str()) {
        return 0.0;
    }

    return v;
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

static double sum_topology_data_yield(const CSV& csv,
                                      const std::vector<std::string>& row,
                                      const std::string& channel,
                                      const std::string& period) {
    double total = 0.0;
    for (const std::string& topo : kTopologies) {
        const int c = col_strict(csv, normalized_raw_yield_col(channel, topo, period, "unpol"));
        total += to_double_or_zero(row[c]);
    }
    return total;
}

static double sum_topology_mc_yield(const CSV& csv,
                                    const std::vector<std::string>& row,
                                    const std::string& channel,
                                    const std::string& period) {
    double total = 0.0;
    for (const std::string& topo : kTopologies) {
        const int c = col_strict(csv, rec_current_corrected_yield_col(channel, topo, period));
        total += to_double_or_zero(row[c]);
    }
    return total;
}

static Triple compute_contamination(double Nmis,
                                    double Npi0_data,
                                    double Npi0_rec_mc,
                                    double Ndvcs_data) {
    Triple out;

    if (!(Nmis > 0.0) || !(Npi0_data > 0.0) || !(Npi0_rec_mc > 0.0) || !(Ndvcs_data > 0.0)) {
        out.v = 0.0;
        out.stat = 0.0;
        out.sys = 0.0;
        return out;
    }

    out.v = (Nmis * Npi0_data) / (Npi0_rec_mc * Ndvcs_data);

    const double rel_var =
        1.0 / Nmis +
        1.0 / Npi0_data +
        1.0 / Npi0_rec_mc +
        1.0 / Ndvcs_data;

    out.stat = std::fabs(out.v) * std::sqrt(std::max(0.0, rel_var));
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

static void compute_period(CSV& csv,
                           const std::vector<RowBin>& bins,
                           const std::string& period,
                           const std::string& output_root_dir) {
    const int c_contam = col_strict(csv, contamination_col(period));

    std::vector<RowContam> points;
    points.resize(csv.rows.size());

    double sum_Ndvcs = 0.0;
    double sum_Nmis = 0.0;
    double sum_Npi0_data = 0.0;
    double sum_Npi0_rec = 0.0;

    for (int r = 0; r < (int)csv.rows.size(); ++r) {
        const RowBin& b = bins[r];

        if (!b.valid) {
            csv.rows[r][c_contam].clear();
            continue;
        }

        const std::vector<std::string>& row = csv.rows[r];

        const double Ndvcs = sum_topology_data_yield(csv, row, "ep->epg", period);
        const double Npi0_data = sum_topology_data_yield(csv, row, "ep->eppi0", period);
        const double Npi0_rec = sum_topology_mc_yield(csv, row, "ep->eppi0", period);
        const double Nmis = sum_topology_mc_yield(csv, row, "ep->eppi0->epg", period);

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

        sum_Ndvcs += Ndvcs;
        sum_Npi0_data += Npi0_data;
        sum_Npi0_rec += Npi0_rec;
        sum_Nmis += Nmis;
    }

    const Triple integrated = compute_contamination(sum_Nmis, sum_Npi0_data, sum_Npi0_rec, sum_Ndvcs);

    std::cout << "[pi0_contamination] " << period
              << " totals: Ndvcs_norm=" << std::setprecision(10) << sum_Ndvcs
              << " Nmis_corr_mc=" << sum_Nmis
              << " Npi0_norm=" << sum_Npi0_data
              << " Npi0_rec_corr_mc=" << sum_Npi0_rec
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
    (void)dvcsDataTrees;
    (void)eppi0DataTrees;
    (void)eppi0RecMcTrees;
    (void)eppi0BkgTrees;
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
