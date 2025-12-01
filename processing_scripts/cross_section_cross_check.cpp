// cross_section_cross_check.cpp
// -----------------------------------------------------------------------------
// Cross-check of experimental cross sections (Hayward pass-2 vs Lee pass-1)
// using *only* CSVs.
// - Lee CSV (pass-1): e.g. imports/all_bin_v3.csv
// - Hayward CSV (pass-2): e.g. output/csvs/dvcs_pass2_analysis.csv
//
// For each valid bin, we compare:
//
//   Lee:     "cross sections, ep->epg, exp"                  (single value)
//   Hayward: "cross sections, ep->epg, exp, Fa18, unpol"    (three-tuple:
//                                                            value, stat_err, syst_err)
//
// We organize the comparison by xB, Q^{2}, and -t. For each (xB, Q^{2}, -t)
// cell, we take all rows matching those ranges and use their provided phiavg
// values as the x-coordinates. We do NOT rebin phi.
//
// Output filenames (per xB index ix):
//   cross_section_counts_xB_<ix>.png
//   cross_section_ratio_xB_<ix>.png
//
// -----------------------------------------------------------------------------

#include "cross_section_cross_check.h"

#include <TCanvas.h>
#include <TGraph.h>
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
#include <cmath>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>
#include <stdexcept>

namespace fs = std::filesystem;

// ---------- small utilities ----------

static inline void info_xs(const std::string& s) {
    std::cout << "[xs] " << s << std::endl;
}

static inline void warn_xs(const std::string& s) {
    std::cout << "[xs][warn] " << s << std::endl;
}

static inline void fatal_xs(const std::string& s) {
    std::cerr << "[xs][FATAL] " << s << std::endl;
    throw std::runtime_error(s);
}

static inline std::string slower_xs(std::string s) {
    for (auto& c : s) {
        c = (char)std::tolower((unsigned char)c);
    }
    return s;
}

// trim leading/trailing whitespace
static inline std::string trim_xs(const std::string& s) {
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

static std::vector<std::string> split_csv_line_xs(const std::string& line) {
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
    out.push_back(cur);
    return out;
}

static std::unordered_map<std::string, int>
build_header_index_xs(const std::vector<std::string>& header) {
    std::unordered_map<std::string, int> m;
    for (int i = 0; i < (int)header.size(); ++i) {
        m[header[i]] = i;
    }
    return m;
}

static const std::string& get_col_ref_xs(const std::vector<std::string>& row,
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

static double ToDouble_xs(const std::string& s) {
    if (s.empty()) return 0.0;
    char* endp = nullptr;
    double v = std::strtod(s.c_str(), &endp);
    if (endp == s.c_str()) return 0.0;
    return v;
}

static int ToInt_xs(const std::string& s) {
    if (s.empty()) return 0;
    char* endp = nullptr;
    long v = std::strtol(s.c_str(), &endp, 10);
    if (endp == s.c_str()) return 0;
    return (int)v;
}

static void require_columns_xs(const std::unordered_map<std::string,int>& idx,
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
        fatal_xs(oss.str());
    }
}

// ---------- Three-tuple parser for Hayward cross section column ----------
// Format: "(value, stat_err, syst_err)" - we extract value and stat_err

struct XsValue_xs {
    double value    = 0.0;
    double stat_err = 0.0;
    bool   valid    = false;
};

static XsValue_xs parse_xs_triplet_xs(const std::string& s) {
    XsValue_xs xv;
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
            parts.push_back(trim_xs(cur));
            cur.clear();
        } else {
            cur.push_back(c);
        }
    }
    parts.push_back(trim_xs(cur));

    if (parts.size() >= 2) {
        xv.value    = ToDouble_xs(parts[0]);
        xv.stat_err = ToDouble_xs(parts[1]);
        xv.valid    = true;
    }

    return xv;
}

// ---------- Lee CSV column detection ----------

struct LeeCsvCols_xs {
    int c_bin;
    int c_xb_min;
    int c_xb_max;
    int c_q2_min;
    int c_q2_max;
    int c_t_min;
    int c_t_max;
    int c_phi_avg;
};

static int find_col_alias_xs(const std::vector<std::string>& header,
                             const std::vector<std::string>& names) {
    for (size_t i = 0; i < header.size(); ++i) {
        std::string h = slower_xs(trim_xs(header[i]));
        for (const auto& raw_name : names) {
            std::string n = slower_xs(trim_xs(raw_name));
            if (!n.empty() && h == n) {
                return (int)i;
            }
        }
        // endfor
    }
    // endfor
    return -1;
}

static LeeCsvCols_xs detect_lee_columns_xs(const std::vector<std::string>& header) {
    if (header.empty()) {
        fatal_xs("Empty header row in Lee CSV");
    }

    LeeCsvCols_xs cols;
    cols.c_bin     = -1;
    cols.c_xb_min  = -1;
    cols.c_xb_max  = -1;
    cols.c_q2_min  = -1;
    cols.c_q2_max  = -1;
    cols.c_t_min   = -1;
    cols.c_t_max   = -1;
    cols.c_phi_avg = -1;

    // 1) Bin index: try named column first, then unlabeled first column
    cols.c_bin = find_col_alias_xs(header, { "bin index", "bin", "idx" });
    if (cols.c_bin < 0) {
        std::string h0 = trim_xs(header[0]);
        if (h0.empty()) {
            cols.c_bin = 0;
        } else {
            std::ostringstream oss;
            oss << "Could not locate bin index column in Lee CSV. "
                << "Tried names: \"bin index\", \"bin\", \"idx\" and unlabeled first column.\n"
                << "Header[0] is \"" << header[0] << "\".";
            fatal_xs(oss.str());
        }
    }

    // 2) xB, Q2, phi
    cols.c_xb_min  = find_col_alias_xs(header, { "xBmin", "xbmin", "xB_min", "xb_min" });
    cols.c_xb_max  = find_col_alias_xs(header, { "xBmax", "xbmax", "xB_max", "xb_max" });
    cols.c_q2_min  = find_col_alias_xs(header, { "Q2min", "q2min", "Q2_min", "q2_min" });
    cols.c_q2_max  = find_col_alias_xs(header, { "Q2max", "q2max", "Q2_max", "q2_max" });
    cols.c_phi_avg = find_col_alias_xs(header, { "phiavg", "phi_avg", "phi_average" });

    // 3) |t| min and max
    cols.c_t_min = find_col_alias_xs(header, { "t_abs_min", "tmin", "t_min" });
    cols.c_t_max = find_col_alias_xs(header, { "t_abs_max", "tmax", "t_max" });

    // 4) Validate required columns
    std::vector<std::string> missing;
    if (cols.c_bin     < 0) missing.push_back("bin index");
    if (cols.c_xb_min  < 0) missing.push_back("xBmin");
    if (cols.c_xb_max  < 0) missing.push_back("xBmax");
    if (cols.c_q2_min  < 0) missing.push_back("Q2min");
    if (cols.c_q2_max  < 0) missing.push_back("Q2max");
    if (cols.c_t_min   < 0) missing.push_back("t_abs_min/tmin");
    if (cols.c_t_max   < 0) missing.push_back("t_abs_max/tmax");
    if (cols.c_phi_avg < 0) missing.push_back("phiavg");

    if (!missing.empty()) {
        std::ostringstream oss;
        oss << "Missing required columns in Lee CSV: ";
        for (size_t i = 0; i < missing.size(); ++i) {
            if (i) oss << ", ";
            oss << "\"" << missing[i] << "\"";
        }
        oss << ".\nHeader row is:";
        for (size_t i = 0; i < header.size(); ++i) {
            oss << "\n  [" << i << "] \"" << header[i] << "\"";
        }
        fatal_xs(oss.str());
    }

    return cols;
}

// ---------- bin / axis structs ----------

struct AxisSets_xs {
    std::vector<std::pair<double,double>> xB;
    std::map<int, std::vector<std::pair<double,double>>> Q2_by_ix;
    std::map<int, std::vector<std::pair<double,double>>> t_by_ix;
};

struct BinRow_xs {
    int    bin_index   = 0;
    double xBmin       = 0.0;
    double xBmax       = 0.0;
    double Q2min       = 0.0;
    double Q2max       = 0.0;
    double tmin        = 0.0;
    double tmax        = 0.0;
    double phiavg      = 0.0;

    double lee_xs      = 0.0;  // Lee cross section (single value, no error)
    double my_xs       = 0.0;  // Hayward cross section (value)
    double my_xs_err   = 0.0;  // Hayward cross section stat error
};

static AxisSets_xs build_axes_from_rows_xs(const std::vector<BinRow_xs>& rows) {
    std::set<std::pair<double,double>> xbset;
    std::map<std::pair<double,double>, std::set<std::pair<double,double>>> q2set_by_xb;
    std::map<std::pair<double,double>, std::set<std::pair<double,double>>> tset_by_xb;

    for (const auto& r : rows) {
        const auto xb = std::make_pair(r.xBmin, r.xBmax);
        xbset.insert(xb);
        q2set_by_xb[xb].insert({r.Q2min, r.Q2max});
        tset_by_xb[xb].insert({r.tmin,  r.tmax});
    }

    AxisSets_xs ax;
    ax.xB.assign(xbset.begin(), xbset.end());
    for (int ix = 0; ix < (int)ax.xB.size(); ++ix) {
        const auto& xb = ax.xB[ix];
        const auto& q2set = q2set_by_xb[xb];
        const auto& tset  = tset_by_xb[xb];
        ax.Q2_by_ix[ix] = { q2set.begin(), q2set.end() };
        ax.t_by_ix[ix]  = { tset.begin(),  tset.end()  };
    }
    return ax;
}

static inline int find_index_xs(const std::pair<double,double>& r,
                                const std::vector<std::pair<double,double>>& v) {
    for (int i = 0; i < (int)v.size(); ++i) {
        if (v[i] == r) return i;
    }
    return -1;
}

// Per-panel data
struct PanelData_xs {
    std::vector<double> phi;
    std::vector<double> val;
    std::vector<double> err;
};

struct PerPanel_xs {
    // key = (ix, iQ, it)
    std::map<std::tuple<int,int,int>, PanelData_xs> lee;
    std::map<std::tuple<int,int,int>, PanelData_xs> hayward;
};

static PerPanel_xs map_to_panels_xs(const std::vector<BinRow_xs>& rows,
                                    const AxisSets_xs& ax) {
    PerPanel_xs pp;

    for (const auto& r : rows) {
        const auto xb = std::make_pair(r.xBmin, r.xBmax);
        const int ix  = find_index_xs(xb, ax.xB);
        if (ix < 0) continue;

        const auto& Q2s = ax.Q2_by_ix.at(ix);
        const auto& Ts  = ax.t_by_ix.at(ix);

        const int iQ = find_index_xs({r.Q2min, r.Q2max}, Q2s);
        const int it = find_index_xs({r.tmin,   r.tmax},  Ts);
        if (iQ < 0 || it < 0) continue;

        auto key = std::make_tuple(ix, iQ, it);

        // Lee data (no error bars)
        if (r.lee_xs > 0.0) {
            pp.lee[key].phi.push_back(r.phiavg);
            pp.lee[key].val.push_back(r.lee_xs);
            pp.lee[key].err.push_back(0.0);
        }

        // Hayward data (with error bars)
        if (r.my_xs > 0.0) {
            pp.hayward[key].phi.push_back(r.phiavg);
            pp.hayward[key].val.push_back(r.my_xs);
            pp.hayward[key].err.push_back(r.my_xs_err);
        }
    }

    // Sort each panel's data by phi
    auto sort_panel = [](PanelData_xs& pd) {
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

    for (auto& kv : pp.lee)     sort_panel(kv.second);
    for (auto& kv : pp.hayward) sort_panel(kv.second);

    return pp;
}

// ---------- plotting helpers ----------

static inline void degreeTicks_xs(double xmin, double ymin, double xmax, double labelSize) {
    TGaxis* ax = new TGaxis(xmin, ymin, xmax, ymin, 0.0, 360.0, 4, "");
    ax->SetLabelFont(42);
    ax->SetLabelSize(labelSize);
    ax->SetLabelOffset(0.012);
    ax->SetTitle("");
    ax->SetTickSize(0.02);
    ax->Draw();
}

static void draw_degree_ticks_here_xs(double labelSize) {
    degreeTicks_xs(gPad->GetUxmin(), gPad->GetUymin(), gPad->GetUxmax(), labelSize);
}

// Markers + error bars
static TGraphErrors* graph_pe1_xs(const std::vector<double>& X,
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

static std::string safe_canvas_name_xs(const std::string& out_png) {
    return fs::path(out_png).filename().string();
}

static void draw_one_canvas_xs(const std::string& title,
                               const std::vector<std::pair<double,double>>& Q2s,
                               const std::vector<std::pair<double,double>>& Ts,
                               const std::function<void(int,int,PanelData_xs&,PanelData_xs&)>& fetchBoth,
                               const std::string& out_png,
                               bool draw_ratio_only) {
    const int nrows = (int)Ts.size();
    const int ncols = (int)Q2s.size();
    if (nrows == 0 || ncols == 0) return;

    const int W = 320 * ncols + 220;
    const int H = 260 * nrows + 260;

    const std::string cname = safe_canvas_name_xs(out_png);
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
        legTop->AddEntry(mRatio, "Hayward/Lee", "p");
        legTop->AddEntry(lnY1,   "y = 1",       "l");
    } else {
        auto* mH = new TMarker(0.0, 0.0, 20);
        mH->SetMarkerColor(kBlack);
        auto* mL = new TMarker(0.0, 0.0, 24);
        mL->SetMarkerColor(kOrange + 7);
        legend_keepalive.push_back(mH);
        legend_keepalive.push_back(mL);
        legTop->AddEntry(mH, "Hayward (pass-2)", "p");
        legTop->AddEntry(mL, "Lee (pass-1)",     "p");
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

            PanelData_xs hayward, lee;
            fetchBoth(ccol, r, hayward, lee);

            double ymin = 0.0;
            double ymax = 1.0;

            if (!draw_ratio_only) {
                // Counts: auto-scale to cover both Lee and Hayward with errors
                ymax = 0.0;
                auto update_ymax = [&](const PanelData_xs& pd) {
                    for (size_t i = 0; i < pd.val.size(); ++i) {
                        double vup = pd.val[i];
                        if (i < pd.err.size()) {
                            vup = pd.val[i] + pd.err[i];
                        }
                        if (vup > ymax) ymax = vup;
                    }
                };
                update_ymax(hayward);
                update_ymax(lee);

                if (ymax <= 0.0) {
                    ymin = 0.0;
                    ymax = 1.0;
                } else {
                    ymin = 0.0;
                    ymax = 1.10 * ymax;
                }

                TH1* frame = gPad->DrawFrame(0.0, ymin, 360.0, ymax);
                frame->GetXaxis()->SetTitle("#phi (deg)");
                frame->GetYaxis()->SetTitle("Cross section");
                frame->GetXaxis()->CenterTitle();
                frame->GetYaxis()->CenterTitle();
                frame->GetXaxis()->SetNdivisions(505);
                frame->GetXaxis()->SetTitleSize(0.060);
                frame->GetYaxis()->SetTitleSize(0.060);
                frame->GetYaxis()->SetLabelSize(0.050);
                frame->GetXaxis()->SetLabelSize(0.050);
                frame->GetXaxis()->SetTitleOffset(1.10);
                frame->GetYaxis()->SetTitleOffset(1.55);

                draw_degree_ticks_here_xs(0.050);

                TLatex lab;
                lab.SetNDC();
                lab.SetTextSize(0.070);
                lab.SetTextAlign(11);
                lab.SetTextFont(42);
                lab.DrawLatex(0.15, 0.92,
                    Form("Q^{2} #in [%.2g, %.2g], -t #in [%.2g, %.2g]",
                         Q2s[ccol].first, Q2s[ccol].second,
                         Ts[r].first,     Ts[r].second));

                // Plot both series
                graph_pe1_xs(hayward.phi, hayward.val, hayward.err, 20, black);  // Hayward with errors
                graph_pe1_xs(lee.phi,     lee.val,     lee.err,     24, orange); // Lee (errors are 0)
            } else {
                // Ratio: compute R = H/L and eR from Hayward error
                const double tol = 20.0;
                std::vector<double> x, y, ey;

                double local_max = 0.0;

                for (size_t i = 0; i < hayward.phi.size(); ++i) {
                    double best_dist = 1e9;
                    int jbest = -1;
                    for (size_t j = 0; j < lee.phi.size(); ++j) {
                        double d = std::fabs(lee.phi[j] - hayward.phi[i]);
                        if (d < best_dist) {
                            best_dist = d;
                            jbest = (int)j;
                        }
                    }
                    if (jbest >= 0 && best_dist <= tol && lee.val[jbest] > 0.0) {
                        double H = hayward.val[i];
                        double L = lee.val[jbest];
                        double R = (H <= 0.0) ? 0.0 : H / L;

                        double sigma_H = (i < hayward.err.size()) ? hayward.err[i] : 0.0;
                        double eR = 0.0;
                        if (H > 0.0 && sigma_H > 0.0) {
                            eR = R * (sigma_H / H);
                        }

                        x.push_back(hayward.phi[i]);
                        y.push_back(R);
                        ey.push_back(eR);

                        if (R + eR > local_max) {
                            local_max = R + eR;
                        }
                    }
                }

                if (local_max <= 0.0) {
                    ymin = 0.0;
                    ymax = 1.0;
                } else {
                    ymin = 0.0;
                    ymax = 1.10 * local_max;
                    if (ymax < 1.0) ymax = 1.0;
                }

                TH1* frame = gPad->DrawFrame(0.0, ymin, 360.0, ymax);
                frame->GetXaxis()->SetTitle("#phi (deg)");
                frame->GetYaxis()->SetTitle("Hayward / Lee");
                frame->GetXaxis()->CenterTitle();
                frame->GetYaxis()->CenterTitle();
                frame->GetXaxis()->SetNdivisions(505);
                frame->GetXaxis()->SetTitleSize(0.060);
                frame->GetYaxis()->SetTitleSize(0.060);
                frame->GetYaxis()->SetLabelSize(0.050);
                frame->GetXaxis()->SetLabelSize(0.050);
                frame->GetXaxis()->SetTitleOffset(1.10);
                frame->GetYaxis()->SetTitleOffset(1.55);

                draw_degree_ticks_here_xs(0.050);

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
                    graph_pe1_xs(x, y, ey, 20, black);

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

// ---------- CSV loaders ----------

static std::vector<BinRow_xs> load_lee_rows_xs(const std::string& lee_csv_path,
                                               std::unordered_map<int,size_t>& bin_to_index) {
    std::ifstream fin(lee_csv_path);
    if (!fin.is_open()) {
        fatal_xs("Cannot open Lee CSV: " + lee_csv_path);
    }

    std::string header_line;
    if (!std::getline(fin, header_line)) {
        fatal_xs("Lee CSV appears empty: " + lee_csv_path);
    }
    std::vector<std::string> header = split_csv_line_xs(header_line);

    LeeCsvCols_xs cols_lee = detect_lee_columns_xs(header);
    auto H = build_header_index_xs(header);

    const std::vector<std::string> required = {
        "valid bin",
        "cross sections, ep->epg, exp"
    };
    require_columns_xs(H, required, "Lee CSV");

    std::vector<BinRow_xs> rows;
    rows.reserve(3000);

    std::string line;
    int input_rows = 0;
    int kept_rows  = 0;

    while (std::getline(fin, line)) {
        ++input_rows;
        if (line.empty()) continue;

        std::vector<std::string> cols = split_csv_line_xs(line);

        const std::string& valid_s = get_col_ref_xs(cols, H, "valid bin");
        int valid = ToInt_xs(valid_s);
        if (valid != 1) continue;

        if ((int)cols.size() <= cols_lee.c_phi_avg) {
            fatal_xs("Row in Lee CSV has fewer columns than expected based on header detection.");
        }

        BinRow_xs r;
        r.bin_index   = ToInt_xs(cols[cols_lee.c_bin]);
        r.xBmin       = ToDouble_xs(cols[cols_lee.c_xb_min]);
        r.xBmax       = ToDouble_xs(cols[cols_lee.c_xb_max]);
        r.Q2min       = ToDouble_xs(cols[cols_lee.c_q2_min]);
        r.Q2max       = ToDouble_xs(cols[cols_lee.c_q2_max]);
        r.tmin        = ToDouble_xs(cols[cols_lee.c_t_min]);
        r.tmax        = ToDouble_xs(cols[cols_lee.c_t_max]);
        r.phiavg      = ToDouble_xs(cols[cols_lee.c_phi_avg]);

        double lee_raw = ToDouble_xs(
            get_col_ref_xs(cols, H, "cross sections, ep->epg, exp")
        );
        r.lee_xs     = lee_raw;
        r.my_xs      = 0.0;
        r.my_xs_err  = 0.0;

        if (bin_to_index.find(r.bin_index) != bin_to_index.end()) {
            fatal_xs("Duplicate bin index in Lee CSV: " + std::to_string(r.bin_index));
        }

        bin_to_index[r.bin_index] = rows.size();
        rows.push_back(r);
        ++kept_rows;
    }

    info_xs("Lee CSV rows read: " + std::to_string(input_rows));
    info_xs("Lee valid rows kept (valid bin == 1): " + std::to_string(kept_rows));
    return rows;
}

static void fill_hayward_xs(const std::string& hayward_csv_path,
                            const std::unordered_map<int,size_t>& bin_to_index,
                            std::vector<BinRow_xs>& rows) {
    std::ifstream fin(hayward_csv_path);
    if (!fin.is_open()) {
        fatal_xs("Cannot open Hayward CSV: " + hayward_csv_path);
    }

    std::string header_line;
    if (!std::getline(fin, header_line)) {
        fatal_xs("Hayward CSV appears empty: " + hayward_csv_path);
    }
    std::vector<std::string> header = split_csv_line_xs(header_line);
    auto H = build_header_index_xs(header);

    const std::vector<std::string> required = {
        "bin index",
        "valid bin",
        "cross sections, ep->epg, exp, Fa18, unpol"
    };
    require_columns_xs(H, required, "Hayward CSV");

    std::string line;
    int input_rows = 0;
    int matched    = 0;

    while (std::getline(fin, line)) {
        ++input_rows;
        if (line.empty()) continue;

        std::vector<std::string> cols = split_csv_line_xs(line);
        const std::string& valid_s = get_col_ref_xs(cols, H, "valid bin");
        int valid = ToInt_xs(valid_s);
        if (valid != 1) continue;

        int bin_index = ToInt_xs(get_col_ref_xs(cols, H, "bin index"));
        auto it = bin_to_index.find(bin_index);
        if (it == bin_to_index.end()) {
            continue;
        }

        XsValue_xs xv = parse_xs_triplet_xs(
            get_col_ref_xs(cols, H, "cross sections, ep->epg, exp, Fa18, unpol")
        );

        BinRow_xs& r = rows[it->second];
        if (xv.valid) {
            r.my_xs     = xv.value;
            r.my_xs_err = xv.stat_err;
        }
        ++matched;
    }

    info_xs("Hayward CSV rows read: " + std::to_string(input_rows));
    info_xs("Hayward valid rows matched to Lee bins: " + std::to_string(matched));
}

// ---------- driver ----------

void plot_cross_section_cross_checks(const std::string& lee_csv_path,
                                     const std::string& hayward_csv_path,
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

    // 1) Load Lee CSV
    std::unordered_map<int,size_t> bin_to_index;
    auto rows = load_lee_rows_xs(lee_csv_path, bin_to_index);

    // 2) Load Hayward CSV
    fill_hayward_xs(hayward_csv_path, bin_to_index, rows);

    // 3) Build axis sets and per-panel maps
    AxisSets_xs ax = build_axes_from_rows_xs(rows);
    PerPanel_xs pp = map_to_panels_xs(rows, ax);

    info_xs("Axis xB bins: " + std::to_string(ax.xB.size()));

    auto make_fetchBoth = [&](int ix) {
        return [&, ix](int iQcol, int irow,
                       PanelData_xs& hayward,
                       PanelData_xs& lee) {
            auto key = std::make_tuple(ix, iQcol, irow);

            auto itH = pp.hayward.find(key);
            if (itH != pp.hayward.end()) {
                hayward = itH->second;
            } else {
                hayward = PanelData_xs();
            }

            auto itL = pp.lee.find(key);
            if (itL != pp.lee.end()) {
                lee = itL->second;
            } else {
                lee = PanelData_xs();
            }
        };
    };

    for (int ix = 0; ix < (int)ax.xB.size(); ++ix) {
        const auto& Q2s = ax.Q2_by_ix[ix];
        const auto& Ts  = ax.t_by_ix[ix];
        if (Q2s.empty() || Ts.empty()) continue;

        const double xb_lo = ax.xB[ix].first;
        const double xb_hi = ax.xB[ix].second;

        const std::string title_counts =
            Form("Cross sections: 10.6 GeV   x_{B} #in [%.3g, %.3g]", xb_lo, xb_hi);
        const std::string title_ratio  =
            Form("Cross section ratio (Hayward/Lee): 10.6 GeV   x_{B} #in [%.3g, %.3g]",
                 xb_lo, xb_hi);

        auto fetchBoth = make_fetchBoth(ix);

        const std::string f_counts =
            (fs::path(output_base_dir) / Form("cross_section_counts_xB_%d.png", ix)).string();
        const std::string f_ratio  =
            (fs::path(output_base_dir) / Form("cross_section_ratio_xB_%d.png",  ix)).string();

        draw_one_canvas_xs(title_counts, Q2s, Ts, fetchBoth, f_counts, /*draw_ratio_only=*/false);
        draw_one_canvas_xs(title_ratio,  Q2s, Ts, fetchBoth, f_ratio,   /*draw_ratio_only=*/true);

        info_xs("Saved: " + f_counts);
        info_xs("Saved: " + f_ratio);
    }
}