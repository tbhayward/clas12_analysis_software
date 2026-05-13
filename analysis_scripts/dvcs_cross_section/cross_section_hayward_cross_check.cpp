// cross_section_hayward_cross_check.cpp
// -----------------------------------------------------------------------------
// Cross-check of experimental cross sections *within* Hayward pass-2,
// comparing compatible run periods / helicities.
//
// Sp18 Inb/Out and Sp18 are unpolarized-only in this workflow. Helicity-
// separated comparisons are restricted to periods with valid helicity-resolved
// luminosities. No separate 10.2 GeV label is used; Sp19 Inb is the 10.2 GeV
// dataset.
//
// Input:
//   - Hayward CSV (pass-2): e.g. output/csvs/dvcs_pass2_analysis.csv
//
// For each valid bin we compare pairs of three-tuples
//   "(value, stat_err, syst_err)"
// in columns like:
//
//   "normed cross sections, ep->epg, exp, Fa18 Inb, unpol"
//   "normed cross sections, ep->epg, exp, Fa18 Out,  unpol"
//   "normed cross sections, ep->epg, exp, Sp18 Inb,  unpol"
//   "normed cross sections, ep->epg, exp, Sp18 Out,  unpol"
//   "normed cross sections, ep->epg, exp, Fa18,      unpol"
//   "normed cross sections, ep->epg, exp, Sp18,      unpol"
//   "normed cross sections, ep->epg, exp, Sp19 Inb,  unpol"
//   "normed cross sections, ep->epg, exp, Fa18 Inb,  pos"
//   "normed cross sections, ep->epg, exp, Fa18 Out,  pos"
//   "normed cross sections, ep->epg, exp, Fa18 Inb,  neg"
//   "normed cross sections, ep->epg, exp, Fa18 Out,  neg"
//
// We organize the comparison by xB, Q^{2}, and -t using the Lee-style
// binning columns in dvcs_pass2_analysis.csv:
//
//   "xBmin", "xBmax", "Q2min", "Q2max", "t_abs_min", "t_abs_max",
//   and a common phi coordinate from:
//   "phiavg, 10.6 GeV"
//
// For each xB bin index ix, we produce:
//
//   cross_section_counts_xB_<ix>.png
//   cross_section_ratio_xB_<ix>.png
//
// in a subdirectory per comparison under output_base_dir.
//
// STANDARD INSTRUCTIONS (DVCS cross-section cross-checks):
//   - For counts canvases (cross sections):
//       * Log y-scale.
//       * Per-xB-bin y-min floors (lowest to highest xB):
//           1e-1, 1e-2, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-4
//       * Per-canvas global y-max (from all subplots, including error bars),
//         used for all subplots on that canvas.
//   - For ratio canvases (A/B):
//       * Linear y-scale, fixed range [0, 3] for all subplots.
// -----------------------------------------------------------------------------

#include "cross_section_hayward_cross_check.h"

#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMarker.h>
#include <TPad.h>
#include <TH1.h>
#include <TStyle.h>
#include <TGaxis.h>
#include <TString.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>
#include <cstdlib>

namespace fs = std::filesystem;

// ---------- small utilities ----------

static inline void info_xs_h(const std::string& s) {
    std::cout << "[xs-hayward] " << s << std::endl;
}

static inline void fatal_xs_h(const std::string& s) {
    std::cerr << "[xs-hayward][FATAL] " << s << std::endl;
    throw std::runtime_error(s);
}

// trim leading/trailing whitespace
static inline std::string trim_xs_h(const std::string& s) {
    size_t i0 = 0;
    while (i0 < s.size() && std::isspace((unsigned char)s[i0])) {
        ++i0;
    }
    size_t i1 = s.size();
    while (i1 > i0 && std::isspace((unsigned char)s[i1 - 1])) {
        --i1;
    }
    return s.substr(i0, i1 - i0);
}

// ---------- CSV helpers ----------

static std::vector<std::string> split_csv_line_xs_h(const std::string& line) {
    std::vector<std::string> out;
    std::string cur;
    cur.reserve(line.size());
    bool in_quotes = false;

    for (char c : line) {
        if (c == '"') {
            in_quotes = !in_quotes;
            continue;
        }
        if (c == ',' && !in_quotes) {
            out.push_back(cur);
            cur.clear();
        } else {
            cur.push_back(c);
        }
    }
    // endfor
    out.push_back(cur);
    return out;
}

static std::unordered_map<std::string,int>
build_header_index_xs_h(const std::vector<std::string>& header) {
    std::unordered_map<std::string,int> m;
    for (int i = 0; i < (int)header.size(); ++i) {
        m[header[i]] = i;
    }
    // endfor
    return m;
}

static const std::string& get_col_ref_xs_h(const std::vector<std::string>& row,
                                           const std::unordered_map<std::string,int>& idx,
                                           const std::string& name) {
    auto it = idx.find(name);
    if (it == idx.end()) {
        static const std::string empty;
        return empty;
    }
    int j = it->second;
    if (j < 0 || j >= (int)row.size()) {
        static const std::string empty;
        return empty;
    }
    return row[j];
}

static double ToDouble_xs_h(const std::string& s) {
    if (s.empty()) return 0.0;
    char* endp = nullptr;
    double v = std::strtod(s.c_str(), &endp);
    if (endp == s.c_str()) return 0.0;
    return v;
}

static int ToInt_xs_h(const std::string& s) {
    if (s.empty()) return 0;
    char* endp = nullptr;
    long v = std::strtol(s.c_str(), &endp, 10);
    if (endp == s.c_str()) return 0;
    return (int)v;
}

static void require_columns_xs_h(const std::unordered_map<std::string,int>& idx,
                                 const std::vector<std::string>& cols,
                                 const std::string& which_csv) {
    std::vector<std::string> missing;
    for (const auto& c : cols) {
        if (idx.find(c) == idx.end()) {
            missing.push_back(c);
        }
        // endfor
    }
    // endfor
    if (!missing.empty()) {
        std::ostringstream oss;
        oss << "Missing columns in " << which_csv << ": ";
        for (size_t i = 0; i < missing.size(); ++i) {
            if (i > 0) oss << ", ";
            oss << '"' << missing[i] << '"';
        }
        fatal_xs_h(oss.str());
    }
}

// ---------- three-tuple parser for Hayward cross section columns ----------
// Format: "(value, stat_err, syst_err)" - we extract value and stat_err

struct XsValue_h {
    double value    = 0.0;
    double stat_err = 0.0;
    bool   valid    = false;
};

static XsValue_h parse_xs_triplet_xs_h(const std::string& s) {
    XsValue_h xv;
    if (s.empty()) return xv;

    size_t p0 = s.find('(');
    size_t p1 = s.find(')');
    if (p0 == std::string::npos || p1 == std::string::npos || p1 <= p0) {
        return xv;
    }

    std::string inner = s.substr(p0 + 1, p1 - p0 - 1);

    std::vector<std::string> parts;
    std::string cur;
    for (char c : inner) {
        if (c == ',') {
            parts.push_back(trim_xs_h(cur));
            cur.clear();
        } else {
            cur.push_back(c);
        }
    }
    // endfor
    parts.push_back(trim_xs_h(cur));

    if (parts.size() >= 2) {
        xv.value    = ToDouble_xs_h(parts[0]);
        xv.stat_err = ToDouble_xs_h(parts[1]);
        xv.valid    = true;
    }

    return xv;
}

// ---------- bin / axis structs ----------

struct BinRow_h {
    int    bin_index   = 0;
    double xBmin       = 0.0;
    double xBmax       = 0.0;
    double Q2min       = 0.0;
    double Q2max       = 0.0;
    double tmin        = 0.0;
    double tmax        = 0.0;
    double phiavg      = 0.0;   // from "phiavg, 10.6 GeV"

    // Cross section values for each column of interest
    std::unordered_map<std::string,XsValue_h> xs;
};

struct AxisSets_h {
    std::vector<std::pair<double,double>> xB;
    std::map<int, std::vector<std::pair<double,double>>> Q2_by_ix;
    std::map<int, std::vector<std::pair<double,double>>> t_by_ix;
};

static AxisSets_h build_axes_from_rows_xs_h(const std::vector<BinRow_h>& rows) {
    std::set<std::pair<double,double>> xbset;
    std::map<std::pair<double,double>, std::set<std::pair<double,double>>> q2set_by_xb;
    std::map<std::pair<double,double>, std::set<std::pair<double,double>>> tset_by_xb;

    for (const auto& r : rows) {
        const auto xb = std::make_pair(r.xBmin, r.xBmax);
        xbset.insert(xb);
        q2set_by_xb[xb].insert({r.Q2min, r.Q2max});
        tset_by_xb[xb].insert({r.tmin,  r.tmax});
    }
    // endfor

    AxisSets_h ax;
    ax.xB.assign(xbset.begin(), xbset.end());
    for (int ix = 0; ix < (int)ax.xB.size(); ++ix) {
        const auto& xb    = ax.xB[ix];
        const auto& q2set = q2set_by_xb[xb];
        const auto& tset  = tset_by_xb[xb];
        ax.Q2_by_ix[ix] = { q2set.begin(), q2set.end() };
        ax.t_by_ix[ix]  = { tset.begin(),  tset.end()  };
    }
    // endfor
    return ax;
}

static inline int find_index_xs_h(const std::pair<double,double>& r,
                                  const std::vector<std::pair<double,double>>& v) {
    for (int i = 0; i < (int)v.size(); ++i) {
        if (v[i] == r) return i;
    }
    // endfor
    return -1;
}

// Per-panel data
struct PanelData_h {
    std::vector<double> phi;
    std::vector<double> val;
    std::vector<double> err;
};

struct PerPanelPair_h {
    // key = (ix, iQ, it)
    std::map<std::tuple<int,int,int>, PanelData_h> A;
    std::map<std::tuple<int,int,int>, PanelData_h> B;
};

// ---------- panel mapper for a given column pair ----------

static PerPanelPair_h map_to_panels_pair_xs_h(const std::vector<BinRow_h>& rows,
                                              const AxisSets_h& ax,
                                              const std::string& colA,
                                              const std::string& colB) {
    PerPanelPair_h pp;

    for (const auto& r : rows) {
        const auto xb = std::make_pair(r.xBmin, r.xBmax);
        const int ix  = find_index_xs_h(xb, ax.xB);
        if (ix < 0) continue;

        const auto& Q2s = ax.Q2_by_ix.at(ix);
        const auto& Ts  = ax.t_by_ix.at(ix);

        const int iQ = find_index_xs_h({r.Q2min, r.Q2max}, Q2s);
        const int it = find_index_xs_h({r.tmin,   r.tmax},  Ts);
        if (iQ < 0 || it < 0) continue;

        auto key = std::make_tuple(ix, iQ, it);

        auto itA = r.xs.find(colA);
        if (itA != r.xs.end() && itA->second.valid && itA->second.value > 0.0) {
            PanelData_h& pdA = pp.A[key];
            pdA.phi.push_back(r.phiavg);
            pdA.val.push_back(itA->second.value);
            pdA.err.push_back(itA->second.stat_err);
        }

        auto itB = r.xs.find(colB);
        if (itB != r.xs.end() && itB->second.valid && itB->second.value > 0.0) {
            PanelData_h& pdB = pp.B[key];
            pdB.phi.push_back(r.phiavg);
            pdB.val.push_back(itB->second.value);
            pdB.err.push_back(itB->second.stat_err);
        }
    }
    // endfor

    // Sort each panel's data by phi
    auto sort_panel = [](PanelData_h& pd) {
        if (pd.phi.empty()) return;
        std::vector<size_t> idx(pd.phi.size());
        for (size_t i = 0; i < idx.size(); ++i) idx[i] = i;
        std::sort(idx.begin(), idx.end(), [&](size_t a, size_t b) {
            return pd.phi[a] < pd.phi[b];
        });
        std::vector<double> phi2, val2, err2;
        phi2.reserve(idx.size());
        val2.reserve(idx.size());
        err2.reserve(idx.size());
        for (size_t i : idx) {
            phi2.push_back(pd.phi[i]);
            val2.push_back(pd.val[i]);
            err2.push_back(pd.err[i]);
        }
        pd.phi = std::move(phi2);
        pd.val = std::move(val2);
        pd.err = std::move(err2);
    };

    for (auto& kv : pp.A) sort_panel(kv.second);
    for (auto& kv : pp.B) sort_panel(kv.second);

    return pp;
}

// ---------- plotting helpers ----------

static inline void degreeTicks_xs_h(double xmin, double ymin, double xmax, double labelSize) {
    TGaxis* ax = new TGaxis(xmin, ymin, xmax, ymin, 0.0, 360.0, 4, "");
    ax->SetLabelFont(42);
    ax->SetLabelSize(labelSize);
    ax->SetLabelOffset(0.012);
    ax->SetTitle("");
    ax->SetTickSize(0.02);
    ax->Draw();
}

static void draw_degree_ticks_here_xs_h(double labelSize) {
    degreeTicks_xs_h(gPad->GetUxmin(), gPad->GetUymin(), gPad->GetUxmax(), labelSize);
}

// Markers + error bars
static TGraphErrors* graph_pe1_xs_h(const std::vector<double>& X,
                                    const std::vector<double>& Y,
                                    const std::vector<double>& EY,
                                    int markerStyle, int color) {
    if (X.empty()) return nullptr;
    std::vector<double> ex(X.size(), 0.0);
    auto* g = new TGraphErrors((int)X.size(),
                               const_cast<double*>(X.data()),
                               const_cast<double*>(Y.data()),
                               ex.data(),
                               const_cast<double*>(EY.data()));
    g->SetMarkerStyle(markerStyle);
    g->SetMarkerColor(color);
    g->SetLineColor(color);
    g->SetLineWidth(1);
    g->Draw("PE1 SAME");
    return g;
}

static std::string safe_canvas_name_xs_h(const std::string& out_png) {
    return fs::path(out_png).filename().string();
}

// OFFICIAL y-min floors for xB bins (lowest to highest xB)
static double get_xb_ymin_floor_xs_h(int ix_xb) {
    static const double floors[] = {
        1e-1,  // xB bin 0
        1e-2,  // xB bin 1
        1e-3,  // xB bin 2
        1e-3,  // xB bin 3
        1e-3,  // xB bin 4
        1e-3,  // xB bin 5
        1e-3,  // xB bin 6
        1e-4   // xB bin 7
    };
    const int n = (int)(sizeof(floors) / sizeof(floors[0]));
    if (ix_xb < 0) {
        return floors[0];
    }
    if (ix_xb >= n) {
        return floors[n - 1];
    }
    return floors[ix_xb];
}

static void draw_one_canvas_pair_xs_h(
    const std::string& title,
    const std::vector<std::pair<double,double>>& Q2s,
    const std::vector<std::pair<double,double>>& Ts,
    const std::function<void(int,int,PanelData_h&,PanelData_h&)>& fetchBoth,
    const std::string& out_png,
    bool draw_ratio_only,
    const std::string& labelA,
    const std::string& labelB,
    int ix_xb
) {
    const int nrows = (int)Ts.size();
    const int ncols = (int)Q2s.size();
    if (nrows == 0 || ncols == 0) return;

    const int W = 320 * ncols + 220;
    const int H = 260 * nrows + 260;

    // -------------------------------------------------------------------------
    // STANDARD INSTRUCTIONS: Precompute global y-range for counts canvases
    // -------------------------------------------------------------------------
    double global_max_counts = 0.0;
    bool   any_positive      = false;
    double y_floor           = 0.0;

    if (!draw_ratio_only) {
        y_floor = get_xb_ymin_floor_xs_h(ix_xb);

        for (int r = 0; r < nrows; ++r) {
            for (int ccol = 0; ccol < ncols; ++ccol) {
                PanelData_h A_tmp, B_tmp;
                fetchBoth(ccol, r, A_tmp, B_tmp);

                auto update_minmax = [&](const PanelData_h& pd) {
                    for (size_t i = 0; i < pd.val.size(); ++i) {
                        double v   = pd.val[i];
                        double es  = (i < pd.err.size()) ? pd.err[i] : 0.0;
                        double vup = v + es;
                        if (vup > global_max_counts) {
                            global_max_counts = vup;
                        }
                        if (vup > 0.0) {
                            any_positive = true;
                        }
                    }
                    // endfor
                };

                update_minmax(A_tmp);
                update_minmax(B_tmp);
            }
        }
    }

    const std::string cname = safe_canvas_name_xs_h(out_png);
    TCanvas* c = new TCanvas(cname.c_str(), cname.c_str(), W, H);
    c->cd();

    // Top band
    TPad* pTop = new TPad("pTop", "pTop", 0.0, 0.90, 1.0, 1.0);
    pTop->SetFillStyle(0);
    pTop->SetBorderSize(0);
    pTop->Draw();
    pTop->cd();

    TLatex head;
    head.SetNDC();
    head.SetTextAlign(22);
    head.SetTextFont(42);
    head.SetTextSize(0.18);
    head.DrawLatex(0.50, 0.65, title.c_str());

    // Legend
    std::vector<TObject*> legend_keepalive;
    TLegend* legTop = new TLegend(0.08, 0.10, 0.92, 0.56);
    legTop->SetNColumns(2);
    legTop->SetBorderSize(0);
    legTop->SetFillStyle(0);
    legTop->SetTextFont(42);
    legTop->SetTextSize(0.22);
    if (draw_ratio_only) {
        auto* mRatio = new TMarker(0.0, 0.0, 20);
        mRatio->SetMarkerColor(kBlack);
        auto* lnY1   = new TLine(0.0, 0.0, 1.0, 0.0);
        lnY1->SetLineStyle(2);
        lnY1->SetLineWidth(2);
        lnY1->SetLineColor(kOrange + 7);
        legend_keepalive.push_back(mRatio);
        legend_keepalive.push_back(lnY1);
        std::string ratioLabel = labelA + "/" + labelB;
        legTop->AddEntry(mRatio, ratioLabel.c_str(), "p");
        legTop->AddEntry(lnY1,   "y = 1",             "l");
    } else {
        auto* mA = new TMarker(0.0, 0.0, 20);
        mA->SetMarkerColor(kBlack);
        auto* mB = new TMarker(0.0, 0.0, 24);
        mB->SetMarkerColor(kOrange + 7);
        legend_keepalive.push_back(mA);
        legend_keepalive.push_back(mB);
        legTop->AddEntry(mA, labelA.c_str(), "p");
        legTop->AddEntry(mB, labelB.c_str(), "p");
    }
    legTop->Draw();

    // Grid
    c->cd();
    TPad* pGrid = new TPad("pGrid", "pGrid", 0.0, 0.00, 1.0, 0.90);
    pGrid->SetFillStyle(0);
    pGrid->SetBorderSize(0);
    pGrid->Draw();
    pGrid->cd();
    pGrid->Divide(ncols, nrows, 0.00, 0.00);

    const int black  = kBlack;
    const int orange = kOrange + 7;

    for (int r = 0; r < nrows; ++r) {
        for (int ccol = 0; ccol < ncols; ++ccol) {
            pGrid->cd(r * ncols + ccol + 1);
            gPad->SetGrid(1, 1);
            gPad->SetTicks(1, 1);
            gPad->SetTopMargin(0.12);
            gPad->SetBottomMargin(0.18);
            gPad->SetLeftMargin(0.18);
            gPad->SetRightMargin(0.08);

            PanelData_h A, B;
            fetchBoth(ccol, r, A, B);

            double ymin = 0.0;
            double ymax = 1.0;

            if (!draw_ratio_only) {
                if (!any_positive || global_max_counts <= 0.0) {
                    gPad->SetLogy(0);
                    ymin = 0.0;
                    ymax = 1.0;
                } else {
                    gPad->SetLogy(1);
                    ymin = y_floor;
                    ymax = 1.10 * global_max_counts;
                    if (ymax <= ymin) {
                        ymax = ymin * 10.0;
                    }
                }

                TH1* frame = gPad->DrawFrame(0.0, ymin, 360.0, ymax);
                frame->GetXaxis()->SetTitle("#phi (deg)");
                frame->GetYaxis()->SetTitle("Normed cross section");
                frame->GetXaxis()->CenterTitle();
                frame->GetYaxis()->CenterTitle();
                frame->GetXaxis()->SetNdivisions(505);
                frame->GetXaxis()->SetTitleSize(0.060);
                frame->GetYaxis()->SetTitleSize(0.060);
                frame->GetYaxis()->SetLabelSize(0.050);
                frame->GetXaxis()->SetLabelSize(0.050);
                frame->GetXaxis()->SetTitleOffset(1.10);
                frame->GetYaxis()->SetTitleOffset(1.55);

                draw_degree_ticks_here_xs_h(0.050);

                TLatex lab;
                lab.SetNDC();
                lab.SetTextSize(0.070);
                lab.SetTextAlign(11);
                lab.SetTextFont(42);
                lab.DrawLatex(0.15, 0.92,
                    Form("Q^{2} #in [%.2g, %.2g], -t #in [%.2g, %.2g]",
                         Q2s[ccol].first, Q2s[ccol].second,
                         Ts[r].first,     Ts[r].second));

                graph_pe1_xs_h(A.phi, A.val, A.err, 20, black);
                graph_pe1_xs_h(B.phi, B.val, B.err, 24, orange);
            } else {
                gPad->SetLogy(0);

                // Compute R = A/B and eR from A stat error only
                const double tol = 20.0;
                std::vector<double> x, y, ey;

                for (size_t i = 0; i < A.phi.size(); ++i) {
                    double best_dist = 1e9;
                    int jbest = -1;
                    for (size_t j = 0; j < B.phi.size(); ++j) {
                        double d = std::fabs(B.phi[j] - A.phi[i]);
                        if (d < best_dist) {
                            best_dist = d;
                            jbest = (int)j;
                        }
                    }
                    // endfor
                    if (jbest >= 0 && best_dist <= tol && B.val[jbest] > 0.0) {
                        double H = A.val[i];
                        double L = B.val[jbest];
                        double R = (H <= 0.0) ? 0.0 : H / L;

                        double sigma_H = (i < A.err.size()) ? A.err[i] : 0.0;
                        double eR = 0.0;
                        if (H > 0.0 && sigma_H > 0.0) {
                            eR = R * (sigma_H / H);
                        }

                        x.push_back(A.phi[i]);
                        y.push_back(R);
                        ey.push_back(eR);
                    }
                }
                // endfor

                ymin = 0.0;
                ymax = 3.0;

                TH1* frame = gPad->DrawFrame(0.0, ymin, 360.0, ymax);
                frame->GetXaxis()->SetTitle("#phi (deg)");
                std::string ytitle = labelA + "/" + labelB;
                frame->GetYaxis()->SetTitle(ytitle.c_str());
                frame->GetXaxis()->CenterTitle();
                frame->GetYaxis()->CenterTitle();
                frame->GetXaxis()->SetNdivisions(505);
                frame->GetXaxis()->SetTitleSize(0.060);
                frame->GetYaxis()->SetTitleSize(0.060);
                frame->GetYaxis()->SetLabelSize(0.050);
                frame->GetXaxis()->SetLabelSize(0.050);
                frame->GetXaxis()->SetTitleOffset(1.10);
                frame->GetYaxis()->SetTitleOffset(1.55);

                draw_degree_ticks_here_xs_h(0.050);

                TLatex lab;
                lab.SetNDC();
                lab.SetTextSize(0.070);
                lab.SetTextAlign(11);
                lab.SetTextFont(42);
                lab.DrawLatex(0.15, 0.92,
                    Form("Q^{2} #in [%.2g, %.2g], -t #in [%.2g, %.2g]",
                         Q2s[ccol].first, Q2s[ccol].second,
                         Ts[r].first,     Ts[r].second));

                if (!x.empty()) {
                    graph_pe1_xs_h(x, y, ey, 20, black);

                    TLine* one = new TLine(0.0, 1.0, 360.0, 1.0);
                    one->SetLineStyle(2);
                    one->SetLineWidth(2);
                    one->SetLineColor(orange);
                    one->Draw("SAME");
                }
            }
        }
    }

    c->SaveAs(out_png.c_str());

    delete legTop;
    delete c;
}

// ---------- CSV loader for Hayward rows ----------

static std::vector<BinRow_h>
load_hayward_rows_xs_h(const std::string& hayward_csv_path,
                       const std::vector<std::string>& xs_columns) {
    std::ifstream fin(hayward_csv_path);
    if (!fin.is_open()) {
        fatal_xs_h("Cannot open Hayward CSV: " + hayward_csv_path);
    }

    std::string header_line;
    if (!std::getline(fin, header_line)) {
        fatal_xs_h("Hayward CSV appears empty: " + hayward_csv_path);
    }
    std::vector<std::string> header = split_csv_line_xs_h(header_line);
    auto H = build_header_index_xs_h(header);

    // Axis + housekeeping columns required
    std::vector<std::string> required = {
        "bin index",
        "valid bin",
        "xBmin",
        "xBmax",
        "Q2min",
        "Q2max",
        "t_abs_min",
        "t_abs_max",
        "phiavg, 10.6 GeV"
    };

    for (const auto& c : xs_columns) {
        required.push_back(c);
    }
    // endfor

    require_columns_xs_h(H, required, "Hayward CSV");

    std::vector<BinRow_h> rows;
    rows.reserve(3000);

    std::string line;
    int input_rows = 0;
    int kept_rows  = 0;

    while (std::getline(fin, line)) {
        ++input_rows;
        if (line.empty()) continue;

        std::vector<std::string> cols = split_csv_line_xs_h(line);

        const std::string& valid_s = get_col_ref_xs_h(cols, H, "valid bin");
        int valid = ToInt_xs_h(valid_s);
        if (valid != 1) continue;

        BinRow_h r;
        r.bin_index = ToInt_xs_h(get_col_ref_xs_h(cols, H, "bin index"));
        r.xBmin     = ToDouble_xs_h(get_col_ref_xs_h(cols, H, "xBmin"));
        r.xBmax     = ToDouble_xs_h(get_col_ref_xs_h(cols, H, "xBmax"));
        r.Q2min     = ToDouble_xs_h(get_col_ref_xs_h(cols, H, "Q2min"));
        r.Q2max     = ToDouble_xs_h(get_col_ref_xs_h(cols, H, "Q2max"));
        r.tmin      = ToDouble_xs_h(get_col_ref_xs_h(cols, H, "t_abs_min"));
        r.tmax      = ToDouble_xs_h(get_col_ref_xs_h(cols, H, "t_abs_max"));
        r.phiavg    = ToDouble_xs_h(get_col_ref_xs_h(cols, H, "phiavg, 10.6 GeV"));

        bool has_any_xs = false;

        for (const auto& cname : xs_columns) {
            const std::string& cell = get_col_ref_xs_h(cols, H, cname);
            XsValue_h xv = parse_xs_triplet_xs_h(cell);
            if (xv.valid) {
                has_any_xs = true;
            }
            r.xs[cname] = xv;
        }
        // endfor

        if (!has_any_xs) continue;

        rows.push_back(r);
        ++kept_rows;
    }
    // end while

    info_xs_h("Hayward CSV rows read: " + std::to_string(input_rows));
    info_xs_h("Hayward valid rows with any cross section: " + std::to_string(kept_rows));
    return rows;
}

// ---------- one comparison driver ----------

static void run_one_comparison_xs_h(
    const std::string& labelA,
    const std::string& labelB,
    const std::string& colA,
    const std::string& colB,
    const std::string& subdir_name,
    const std::vector<BinRow_h>& rows,
    const AxisSets_h& ax,
    const std::string& output_base_dir
) {
    fs::path out_dir = fs::path(output_base_dir) / subdir_name;
    fs::create_directories(out_dir);

    info_xs_h("Running cross-section comparison: " + labelA + " vs " + labelB +
              " -> " + out_dir.string());

    PerPanelPair_h pp = map_to_panels_pair_xs_h(rows, ax, colA, colB);

    auto make_fetchBoth = [&](int ix) {
        return [&, ix](int iQcol, int irow,
                       PanelData_h& A,
                       PanelData_h& B) {
            auto key = std::make_tuple(ix, iQcol, irow);

            auto itA = pp.A.find(key);
            if (itA != pp.A.end()) {
                A = itA->second;
            } else {
                A = PanelData_h();
            }

            auto itB = pp.B.find(key);
            if (itB != pp.B.end()) {
                B = itB->second;
            } else {
                B = PanelData_h();
            }
        };
    };

    for (int ix = 0; ix < (int)ax.xB.size(); ++ix) {
        const auto& Q2s = ax.Q2_by_ix.at(ix);
        const auto& Ts  = ax.t_by_ix.at(ix);
        if (Q2s.empty() || Ts.empty()) continue;

        const double xb_lo = ax.xB[ix].first;
        const double xb_hi = ax.xB[ix].second;

        std::string title_counts = Form(
            "Normed cross sections: 10.6 GeV   %s vs %s   x_{B} #in [%.3g, %.3g]",
            labelA.c_str(), labelB.c_str(), xb_lo, xb_hi
        );
        std::string title_ratio = Form(
            "Normed cross section ratio: 10.6 GeV   %s/%s   x_{B} #in [%.3g, %.3g]",
            labelA.c_str(), labelB.c_str(), xb_lo, xb_hi
        );

        auto fetchBoth = make_fetchBoth(ix);

        fs::path f_counts = out_dir / Form("cross_section_counts_xB_%d.png", ix);
        fs::path f_ratio  = out_dir / Form("cross_section_ratio_xB_%d.png",  ix);

        draw_one_canvas_pair_xs_h(
            title_counts, Q2s, Ts, fetchBoth, f_counts.string(),
            /*draw_ratio_only=*/false, labelA, labelB, ix
        );
        draw_one_canvas_pair_xs_h(
            title_ratio,  Q2s, Ts, fetchBoth, f_ratio.string(),
            /*draw_ratio_only=*/true,  labelA, labelB, ix
        );

        info_xs_h("Saved: " + f_counts.string());
        info_xs_h("Saved: " + f_ratio.string());
    }
    // endfor
}

// ---------- top-level driver ----------

void plot_cross_section_hayward_cross_checks(const std::string& hayward_csv_path,
                                             const std::string& output_base_dir) {
    fs::create_directories(output_base_dir);

    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetLineWidth(2);
    gStyle->SetFrameLineWidth(2);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetLabelFont(42, "XYZ");
    gStyle->SetTextFont(42);

    // All unique *normed* cross section columns we will require
    const std::vector<std::string> xs_columns = {
        // unpolarized
        "normed cross sections, ep->epg, exp, Fa18 Inb, unpol",
        "normed cross sections, ep->epg, exp, Fa18 Out, unpol",
        "normed cross sections, ep->epg, exp, Sp18 Inb, unpol",
        "normed cross sections, ep->epg, exp, Sp18 Out, unpol",
        "normed cross sections, ep->epg, exp, Fa18, unpol",
        "normed cross sections, ep->epg, exp, Sp18, unpol",
        "normed cross sections, ep->epg, exp, Sp19 Inb, unpol",
        // helicity-separated Fa18
        "normed cross sections, ep->epg, exp, Fa18 Inb, pos",
        "normed cross sections, ep->epg, exp, Fa18 Out, pos",
        "normed cross sections, ep->epg, exp, Fa18 Inb, neg",
        "normed cross sections, ep->epg, exp, Fa18 Out, neg"
    };

    // 1) Load Hayward rows (axis info + all needed normed cross-section columns)
    auto rows = load_hayward_rows_xs_h(hayward_csv_path, xs_columns);

    // 2) Build axis sets from these rows
    AxisSets_h ax = build_axes_from_rows_xs_h(rows);

    info_xs_h("Axis xB bins: " + std::to_string(ax.xB.size()));

    struct Comparison {
        std::string labelA;
        std::string labelB;
        std::string colA;
        std::string colB;
        std::string subdir;
    };

    std::vector<Comparison> comps;

    // Fa18 Inb vs Fa18 Out (unpolarized)
    comps.push_back({
        "Fa18 Inb (unpol)",
        "Fa18 Out (unpol)",
        "normed cross sections, ep->epg, exp, Fa18 Inb, unpol",
        "normed cross sections, ep->epg, exp, Fa18 Out, unpol",
        "Fa18Inb_vs_Fa18Out_unpol"
    });

    // Fa18 Inb vs Sp18 Inb (unpolarized)
    comps.push_back({
        "Fa18 Inb (unpol)",
        "Sp18 Inb (unpol)",
        "normed cross sections, ep->epg, exp, Fa18 Inb, unpol",
        "normed cross sections, ep->epg, exp, Sp18 Inb, unpol",
        "Fa18Inb_vs_Sp18Inb_unpol"
    });

    // Fa18 Out vs Sp18 Out (unpolarized)
    comps.push_back({
        "Fa18 Out (unpol)",
        "Sp18 Out (unpol)",
        "normed cross sections, ep->epg, exp, Fa18 Out, unpol",
        "normed cross sections, ep->epg, exp, Sp18 Out, unpol",
        "Fa18Out_vs_Sp18Out_unpol"
    });

    // Sp18 Inb vs Sp18 Out (unpolarized)
    comps.push_back({
        "Sp18 Inb (unpol)",
        "Sp18 Out (unpol)",
        "normed cross sections, ep->epg, exp, Sp18 Inb, unpol",
        "normed cross sections, ep->epg, exp, Sp18 Out, unpol",
        "Sp18Inb_vs_Sp18Out_unpol"
    });

    // Fa18 vs Sp18 (unpolarized combined)
    comps.push_back({
        "Fa18 (unpol)",
        "Sp18 (unpol)",
        "normed cross sections, ep->epg, exp, Fa18, unpol",
        "normed cross sections, ep->epg, exp, Sp18, unpol",
        "Fa18_vs_Sp18_unpol"
    });

    // Fa18 Inb vs Sp19 Inb (unpolarized)
    comps.push_back({
        "Fa18 Inb (unpol)",
        "Sp19 Inb (unpol)",
        "normed cross sections, ep->epg, exp, Fa18 Inb, unpol",
        "normed cross sections, ep->epg, exp, Sp19 Inb, unpol",
        "Fa18Inb_vs_Sp19Inb_unpol"
    });

    // Fa18 Inb vs Fa18 Out (pos helicity)
    comps.push_back({
        "Fa18 Inb (pos)",
        "Fa18 Out (pos)",
        "normed cross sections, ep->epg, exp, Fa18 Inb, pos",
        "normed cross sections, ep->epg, exp, Fa18 Out, pos",
        "Fa18Inb_vs_Fa18Out_pos"
    });

    // Fa18 Inb vs Fa18 Out (neg helicity)
    comps.push_back({
        "Fa18 Inb (neg)",
        "Fa18 Out (neg)",
        "normed cross sections, ep->epg, exp, Fa18 Inb, neg",
        "normed cross sections, ep->epg, exp, Fa18 Out, neg",
        "Fa18Inb_vs_Fa18Out_neg"
    });

    for (const auto& cmp : comps) {
        run_one_comparison_xs_h(
            cmp.labelA, cmp.labelB,
            cmp.colA,   cmp.colB,
            cmp.subdir,
            rows, ax,
            output_base_dir
        );
    }
}