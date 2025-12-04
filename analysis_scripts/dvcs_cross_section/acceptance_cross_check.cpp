// acceptance_cross_check.cpp
// -----------------------------------------------------------------------------
// Cross-check of acceptances (Hayward pass-2 vs Lee pass-1) using *only* CSVs.
// - Lee CSV (pass-1):  imports/all_bin_v3.csv
// - Hayward CSV (pass-2): output/csvs/dvcs_pass2_analysis.csv
//
// For each valid bin, we compare:
//
//   Lee:     "acceptance, ep->epg, sim(KM15)"     (single value)
//   Hayward: "acceptance, Fa18 Inb"               (three-tuple: value, stat_err, syst_err)
//             "acceptance, Fa18 Out"              (three-tuple: value, stat_err, syst_err)
//
// We organize the comparison by xB, Q^{2}, and -t. For each (xB, Q^{2}, -t)
// cell, we take all rows matching those ranges and use their provided phiavg
// values as the x-coordinates. We do NOT rebin phi.
//
// For each xB bin we produce two canvases:
//   - acceptance_counts_xB_<ix>.png  : Lee + Fa18 Inb + Fa18 Out
//   - acceptance_ratio_xB_<ix>.png   : (Fa18 Inb / Lee) and (Fa18 Out / Lee)
//
// Within each canvas, all subplots share the same y-axis range for consistency.
// -----------------------------------------------------------------------------

#include "acceptance_cross_check.h"

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

static inline void info_acc(const std::string& s) {
    std::cout << "[acc] " << s << std::endl;
}

static inline void warn_acc(const std::string& s) {
    std::cout << "[acc][warn] " << s << std::endl;
}

static inline void fatal_acc(const std::string& s) {
    std::cerr << "[acc][FATAL] " << s << std::endl;
    throw std::runtime_error(s);
}

static inline std::string slower_acc(std::string s) {
    for (auto& c : s) {
        c = (char)std::tolower((unsigned char)c);
    }
    return s;
}

// trim leading/trailing whitespace
static inline std::string trim_acc(const std::string& s) {
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

static std::vector<std::string> split_csv_line_acc(const std::string& line) {
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
build_header_index_acc(const std::vector<std::string>& header) {
    std::unordered_map<std::string, int> m;
    for (int i = 0; i < (int)header.size(); ++i) {
        m[header[i]] = i;
    }
    return m;
}

static const std::string& get_col_ref_acc(const std::vector<std::string>& row,
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

static double ToDouble_acc(const std::string& s) {
    if (s.empty()) return 0.0;
    char* endp = nullptr;
    double v = std::strtod(s.c_str(), &endp);
    if (endp == s.c_str()) return 0.0;
    return v;
}

static int ToInt_acc(const std::string& s) {
    if (s.empty()) return 0;
    char* endp = nullptr;
    long v = std::strtol(s.c_str(), &endp, 10);
    if (endp == s.c_str()) return 0;
    return (int)v;
}

static void require_columns_acc(const std::unordered_map<std::string,int>& idx,
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
        fatal_acc(oss.str());
    }
}

// ---------- Three-tuple parser for Hayward acceptance columns ----------
// Format: "(value, stat_err, syst_err)"
// We extract:
//   - value    = parts[0]
//   - stat_err = parts[1]

struct AccValue_acc {
    double value    = 0.0;
    double stat_err = 0.0;
    bool   valid    = false;
};

static AccValue_acc parse_acc_triplet_acc(const std::string& s) {
    AccValue_acc av;
    if (s.empty()) return av;

    size_t p0 = s.find('(');
    size_t p1 = s.find(')');
    if (p0 == std::string::npos || p1 == std::string::npos || p1 <= p0) {
        return av;
    }

    std::string inner = s.substr(p0 + 1, p1 - p0 - 1);

    std::vector<std::string> parts;
    std::string cur;
    for (char c : inner) {
        if (c == ',') {
            parts.push_back(trim_acc(cur));
            cur.clear();
        } else {
            cur.push_back(c);
        }
    }
    parts.push_back(trim_acc(cur));

    if (parts.size() >= 2) {
        av.value    = ToDouble_acc(parts[0]);
        av.stat_err = ToDouble_acc(parts[1]);
        av.valid    = true;
    }

    return av;
}

// ---------- Lee CSV column detection ----------

struct LeeCsvCols_acc {
    int c_bin;
    int c_xb_min;
    int c_xb_max;
    int c_q2_min;
    int c_q2_max;
    int c_t_min;
    int c_t_max;
    int c_phi_avg;
};

static int find_col_alias_acc(const std::vector<std::string>& header,
                              const std::vector<std::string>& names) {
    for (size_t i = 0; i < header.size(); ++i) {
        std::string h = slower_acc(trim_acc(header[i]));
        for (const auto& raw_name : names) {
            std::string n = slower_acc(trim_acc(raw_name));
            if (!n.empty() && h == n) {
                return (int)i;
            }
        }
        // endfor
    }
    // endfor
    return -1;
}

static LeeCsvCols_acc detect_lee_columns_acc(const std::vector<std::string>& header) {
    if (header.empty()) {
        fatal_acc("Empty header row in Lee CSV");
    }

    LeeCsvCols_acc cols;
    cols.c_bin     = -1;
    cols.c_xb_min  = -1;
    cols.c_xb_max  = -1;
    cols.c_q2_min  = -1;
    cols.c_q2_max  = -1;
    cols.c_t_min   = -1;
    cols.c_t_max   = -1;
    cols.c_phi_avg = -1;

    // 1) Bin index: try named column first, then unlabeled first column
    cols.c_bin = find_col_alias_acc(header, { "bin index", "bin", "idx" });
    if (cols.c_bin < 0) {
        std::string h0 = trim_acc(header[0]);
        if (h0.empty()) {
            cols.c_bin = 0;
        } else {
            std::ostringstream oss;
            oss << "Could not locate bin index column in Lee CSV. "
                << "Tried names: \"bin index\", \"bin\", \"idx\" and unlabeled first column.\n"
                << "Header[0] is \"" << header[0] << "\".";
            fatal_acc(oss.str());
        }
    }

    // 2) xB, Q2, phi
    cols.c_xb_min  = find_col_alias_acc(header, { "xBmin", "xbmin", "xB_min", "xb_min" });
    cols.c_xb_max  = find_col_alias_acc(header, { "xBmax", "xbmax", "xB_max", "xb_max" });
    cols.c_q2_min  = find_col_alias_acc(header, { "Q2min", "q2min", "Q2_min", "q2_min" });
    cols.c_q2_max  = find_col_alias_acc(header, { "Q2max", "q2max", "Q2_max", "q2_max" });
    cols.c_phi_avg = find_col_alias_acc(header, { "phiavg", "phi_avg", "phi_average" });

    // 3) |t| min and max
    cols.c_t_min = find_col_alias_acc(header, { "t_abs_min", "tmin", "t_min" });
    cols.c_t_max = find_col_alias_acc(header, { "t_abs_max", "tmax", "t_max" });

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
        fatal_acc(oss.str());
    }

    return cols;
}

// ---------- bin / axis structs ----------

struct AxisSets_acc {
    std::vector<std::pair<double,double>> xB;
    std::map<int, std::vector<std::pair<double,double>>> Q2_by_ix;
    std::map<int, std::vector<std::pair<double,double>>> t_by_ix;
};

struct BinRow_acc {
    int    bin_index       = 0;
    double xBmin           = 0.0;
    double xBmax           = 0.0;
    double Q2min           = 0.0;
    double Q2max           = 0.0;
    double tmin            = 0.0;
    double tmax            = 0.0;
    double phiavg          = 0.0;

    double lee_acc         = 0.0;  // Lee acceptance (single value, no error)
    double my_acc_inb      = 0.0;  // Hayward Fa18 Inb acceptance
    double my_acc_inb_err  = 0.0;  // Hayward Fa18 Inb stat error
    double my_acc_out      = 0.0;  // Hayward Fa18 Out acceptance
    double my_acc_out_err  = 0.0;  // Hayward Fa18 Out stat error
};

static AxisSets_acc build_axes_from_rows_acc(const std::vector<BinRow_acc>& rows) {
    std::set<std::pair<double,double>> xbset;
    std::map<std::pair<double,double>, std::set<std::pair<double,double>>> q2set_by_xb;
    std::map<std::pair<double,double>, std::set<std::pair<double,double>>> tset_by_xb;

    for (const auto& r : rows) {
        const auto xb = std::make_pair(r.xBmin, r.xBmax);
        xbset.insert(xb);
        q2set_by_xb[xb].insert({r.Q2min, r.Q2max});
        tset_by_xb[xb].insert({r.tmin,  r.tmax});
    }

    AxisSets_acc ax;
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

static inline int find_index_acc(const std::pair<double,double>& r,
                                 const std::vector<std::pair<double,double>>& v) {
    for (int i = 0; i < (int)v.size(); ++i) {
        if (v[i] == r) return i;
    }
    return -1;
}

// Per-panel data
struct PanelData_acc {
    std::vector<double> phi;
    std::vector<double> val;
    std::vector<double> err;
};

struct PerPanel_acc {
    // key = (ix, iQ, it)
    std::map<std::tuple<int,int,int>, PanelData_acc> lee;
    std::map<std::tuple<int,int,int>, PanelData_acc> inb;
    std::map<std::tuple<int,int,int>, PanelData_acc> out;
};

static PerPanel_acc map_to_panels_acc(const std::vector<BinRow_acc>& rows,
                                      const AxisSets_acc& ax) {
    PerPanel_acc pp;

    for (const auto& r : rows) {
        const auto xb = std::make_pair(r.xBmin, r.xBmax);
        const int ix  = find_index_acc(xb, ax.xB);
        if (ix < 0) continue;

        const auto& Q2s = ax.Q2_by_ix.at(ix);
        const auto& Ts  = ax.t_by_ix.at(ix);

        const int iQ = find_index_acc({r.Q2min, r.Q2max}, Q2s);
        const int it = find_index_acc({r.tmin,   r.tmax},  Ts);
        if (iQ < 0 || it < 0) continue;

        auto key = std::make_tuple(ix, iQ, it);

        // Lee data (no error bars)
        if (r.lee_acc > 0.0) {
            pp.lee[key].phi.push_back(r.phiavg);
            pp.lee[key].val.push_back(r.lee_acc);
            pp.lee[key].err.push_back(0.0);
        }

        // Hayward Fa18 Inb
        if (r.my_acc_inb > 0.0) {
            pp.inb[key].phi.push_back(r.phiavg);
            pp.inb[key].val.push_back(r.my_acc_inb);
            pp.inb[key].err.push_back(r.my_acc_inb_err);
        }

        // Hayward Fa18 Out
        if (r.my_acc_out > 0.0) {
            pp.out[key].phi.push_back(r.phiavg);
            pp.out[key].val.push_back(r.my_acc_out);
            pp.out[key].err.push_back(r.my_acc_out_err);
        }
    }

    // Sort each panel's data by phi
    auto sort_panel = [](PanelData_acc& pd) {
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

    for (auto& kv : pp.lee) sort_panel(kv.second);
    for (auto& kv : pp.inb) sort_panel(kv.second);
    for (auto& kv : pp.out) sort_panel(kv.second);

    return pp;
}

// ---------- plotting helpers ----------

static inline void degreeTicks_acc(double xmin, double ymin, double xmax, double labelSize) {
    TGaxis* ax = new TGaxis(xmin, ymin, xmax, ymin, 0.0, 360.0, 4, "");
    ax->SetLabelFont(42);
    ax->SetLabelSize(labelSize);
    ax->SetLabelOffset(0.012);
    ax->SetTitle("");
    ax->SetTickSize(0.02);
    ax->Draw();
}

static void draw_degree_ticks_here_acc(double labelSize) {
    degreeTicks_acc(gPad->GetUxmin(), gPad->GetUymin(), gPad->GetUxmax(), labelSize);
}

// Markers + error bars
static TGraphErrors* graph_pe1_acc(const std::vector<double>& X,
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

static std::string safe_canvas_name_acc(const std::string& out_png) {
    return fs::path(out_png).filename().string();
}

// Helper to compute the maximum y (including error) for counts in one panel
static double panel_counts_max_acc(const PanelData_acc& lee,
                                   const PanelData_acc& inb,
                                   const PanelData_acc& out) {
    double ymax = 0.0;
    auto update = [&](const PanelData_acc& pd) {
        for (size_t i = 0; i < pd.val.size(); ++i) {
            double vup = pd.val[i];
            if (i < pd.err.size()) {
                vup = pd.val[i] + pd.err[i];
            }
            if (vup > ymax) ymax = vup;
        }
    };
    update(lee);
    update(inb);
    update(out);
    return ymax;
}

// Helper to compute max ratio (including error) for a single Hayward series vs Lee
static double panel_ratio_max_one_acc(const PanelData_acc& pd,
                                      const PanelData_acc& lee,
                                      double tol_deg) {
    double local_max = 0.0;
    for (size_t i = 0; i < pd.phi.size(); ++i) {
        double best_dist = 1e9;
        int jbest = -1;
        for (size_t j = 0; j < lee.phi.size(); ++j) {
            double d = std::fabs(lee.phi[j] - pd.phi[i]);
            if (d < best_dist) {
                best_dist = d;
                jbest = (int)j;
            }
        }
        if (jbest >= 0 && best_dist <= tol_deg && lee.val[jbest] > 0.0) {
            double H = pd.val[i];
            double L = lee.val[jbest];
            double R = (H <= 0.0) ? 0.0 : H / L;

            double sigma_H = (i < pd.err.size()) ? pd.err[i] : 0.0;
            double eR = 0.0;
            if (H > 0.0 && sigma_H > 0.0) {
                eR = R * (sigma_H / H);
            }

            if (R + eR > local_max) {
                local_max = R + eR;
            }
        }
    }
    return local_max;
}

// Main canvas drawing function with per-canvas y-range standardization
static void draw_one_canvas_acc(
    const std::string& title,
    const std::vector<std::pair<double,double>>& Q2s,
    const std::vector<std::pair<double,double>>& Ts,
    const std::function<void(int,int,PanelData_acc&,PanelData_acc&,PanelData_acc&)>& fetchAll,
    const std::string& out_png,
    bool draw_ratio_only
) {
    const int nrows = (int)Ts.size();
    const int ncols = (int)Q2s.size();
    if (nrows == 0 || ncols == 0) return;

    const int W = 320 * ncols + 220;
    const int H = 260 * nrows + 260;

    const std::string cname = safe_canvas_name_acc(out_png);
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

    const int color_inb  = kBlack;
    const int color_out  = kBlue + 1;
    const int color_lee  = kOrange + 7;

    if (draw_ratio_only) {
        auto* mInb = new TMarker(0.0, 0.0, 20);
        mInb->SetMarkerColor(color_inb);
        auto* mOut = new TMarker(0.0, 0.0, 21);
        mOut->SetMarkerColor(color_out);
        auto* lnY1 = new TLine(0.0, 0.0, 1.0, 0.0);
        lnY1->SetLineStyle(2);
        lnY1->SetLineWidth(2);
        lnY1->SetLineColor(color_lee);

        legend_keepalive.push_back(mInb);
        legend_keepalive.push_back(mOut);
        legend_keepalive.push_back(lnY1);

        legTop->AddEntry(mInb, "Fa18 Inb / Lee",  "p");
        legTop->AddEntry(mOut, "Fa18 Out / Lee",  "p");
        legTop->AddEntry(lnY1, "y = 1",           "l");
    } else {
        auto* mInb = new TMarker(0.0, 0.0, 20);
        mInb->SetMarkerColor(color_inb);
        auto* mOut = new TMarker(0.0, 0.0, 21);
        mOut->SetMarkerColor(color_out);
        auto* mLee = new TMarker(0.0, 0.0, 24);
        mLee->SetMarkerColor(color_lee);

        legend_keepalive.push_back(mInb);
        legend_keepalive.push_back(mOut);
        legend_keepalive.push_back(mLee);

        legTop->AddEntry(mInb, "Hayward Fa18 Inb", "p");
        legTop->AddEntry(mOut, "Hayward Fa18 Out", "p");
        legTop->AddEntry(mLee, "Lee (pass-1)",     "p");
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

    // ---------- First pass: determine per-canvas y-range ----------
    double global_ymin = 0.0;
    double global_ymax = 0.0;
    const double phi_match_tol = 20.0; // degrees, for ratio matching

    if (!draw_ratio_only) {
        // Counts: scan all panels for max(value + err)
        double max_counts = 0.0;
        for (int r = 0; r < nrows; ++r) {
            for (int ccol = 0; ccol < ncols; ++ccol) {
                PanelData_acc inb, out, lee;
                fetchAll(ccol, r, inb, out, lee);
                double local_max = panel_counts_max_acc(lee, inb, out);
                if (local_max > max_counts) max_counts = local_max;
            }
        }
        if (max_counts <= 0.0) {
            global_ymin = 0.0;
            global_ymax = 1.0;
        } else {
            global_ymin = 0.0;
            global_ymax = 1.10 * max_counts;
        }
    } else {
        // Ratios: scan all panels for max(R + eR) for Inb and Out
        double max_ratio = 0.0;
        for (int r = 0; r < nrows; ++r) {
            for (int ccol = 0; ccol < ncols; ++ccol) {
                PanelData_acc inb, out, lee;
                fetchAll(ccol, r, inb, out, lee);

                double local_inb = panel_ratio_max_one_acc(inb, lee, phi_match_tol);
                double local_out = panel_ratio_max_one_acc(out, lee, phi_match_tol);

                if (local_inb > max_ratio) max_ratio = local_inb;
                if (local_out > max_ratio) max_ratio = local_out;
            }
        }
        if (max_ratio <= 0.0) {
            global_ymin = 0.0;
            global_ymax = 1.0;
        } else {
            global_ymin = 0.0;
            global_ymax = 1.10 * max_ratio;
            if (global_ymax < 1.0) {
                global_ymax = 1.0;
            }
        }
    }

    // ---------- Second pass: draw all pads with common y-range ----------
    for (int r = 0; r < nrows; ++r) {
        for (int ccol = 0; ccol < ncols; ++ccol) {
            pGrid->cd(r * ncols + ccol + 1);
            gPad->SetGrid(1, 1);
            gPad->SetTicks(1, 1);
            gPad->SetTopMargin(0.12);
            gPad->SetBottomMargin(0.18);
            gPad->SetLeftMargin(0.18);
            gPad->SetRightMargin(0.08);

            PanelData_acc inb, out, lee;
            fetchAll(ccol, r, inb, out, lee);

            if (!draw_ratio_only) {
                TH1* frame = gPad->DrawFrame(0.0, global_ymin, 360.0, global_ymax);
                frame->GetXaxis()->SetTitle("#phi (deg)");
                frame->GetYaxis()->SetTitle("Acceptance");
                frame->GetXaxis()->CenterTitle();
                frame->GetYaxis()->CenterTitle();
                frame->GetXaxis()->SetNdivisions(505);
                frame->GetXaxis()->SetTitleSize(0.060);
                frame->GetYaxis()->SetTitleSize(0.060);
                frame->GetYaxis()->SetLabelSize(0.050);
                frame->GetXaxis()->SetLabelSize(0.050);
                frame->GetXaxis()->SetTitleOffset(1.10);
                frame->GetYaxis()->SetTitleOffset(1.55);

                draw_degree_ticks_here_acc(0.050);

                TLatex lab;
                lab.SetNDC();
                lab.SetTextSize(0.070);
                lab.SetTextAlign(11);
                lab.SetTextFont(42);
                lab.DrawLatex(0.15, 0.92,
                    Form("Q^{2} #in [%.2g, %.2g], -t #in [%.2g, %.2g]",
                         Q2s[ccol].first, Q2s[ccol].second,
                         Ts[r].first,     Ts[r].second));

                // Plot all three series
                graph_pe1_acc(inb.phi, inb.val, inb.err, 20, color_inb);
                graph_pe1_acc(out.phi, out.val, out.err, 21, color_out);
                graph_pe1_acc(lee.phi, lee.val, lee.err, 24, color_lee);
            } else {
                // Ratio: compute arrays for Inb/Lee and Out/Lee and draw in common range
                std::vector<double> x_inb, y_inb, ey_inb;
                std::vector<double> x_out, y_out, ey_out;

                auto fill_ratio = [&](const PanelData_acc& pd,
                                      std::vector<double>& X,
                                      std::vector<double>& Y,
                                      std::vector<double>& EY) {
                    for (size_t i = 0; i < pd.phi.size(); ++i) {
                        double best_dist = 1e9;
                        int jbest = -1;
                        for (size_t j = 0; j < lee.phi.size(); ++j) {
                            double d = std::fabs(lee.phi[j] - pd.phi[i]);
                            if (d < best_dist) {
                                best_dist = d;
                                jbest = (int)j;
                            }
                        }
                        if (jbest >= 0 && best_dist <= phi_match_tol && lee.val[jbest] > 0.0) {
                            double H = pd.val[i];
                            double L = lee.val[jbest];
                            double R = (H <= 0.0) ? 0.0 : H / L;

                            double sigma_H = (i < pd.err.size()) ? pd.err[i] : 0.0;
                            double eR = 0.0;
                            if (H > 0.0 && sigma_H > 0.0) {
                                eR = R * (sigma_H / H);
                            }

                            X.push_back(pd.phi[i]);
                            Y.push_back(R);
                            EY.push_back(eR);
                        }
                    }
                };

                fill_ratio(inb, x_inb, y_inb, ey_inb);
                fill_ratio(out, x_out, y_out, ey_out);

                TH1* frame = gPad->DrawFrame(0.0, global_ymin, 360.0, global_ymax);
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

                draw_degree_ticks_here_acc(0.050);

                TLatex lab;
                lab.SetNDC();
                lab.SetTextSize(0.070);
                lab.SetTextAlign(11);
                lab.SetTextFont(42);
                lab.DrawLatex(0.15, 0.92,
                    Form("Q^{2} #in [%.2g, %.2g], -t #in [%.2g, %.2g]",
                         Q2s[ccol].first, Q2s[ccol].second,
                         Ts[r].first,     Ts[r].second));

                if (!x_inb.empty()) {
                    graph_pe1_acc(x_inb, y_inb, ey_inb, 20, color_inb);
                }
                if (!x_out.empty()) {
                    graph_pe1_acc(x_out, y_out, ey_out, 21, color_out);
                }

                TLine* one = new TLine(0.0, 1.0, 360.0, 1.0);
                one->SetLineStyle(2);
                one->SetLineWidth(2);
                one->SetLineColor(color_lee);
                one->Draw("SAME");
            }
        }
    }

    c->SaveAs(out_png.c_str());

    delete legTop;
    delete c;
}

// ---------- CSV loaders ----------

static std::vector<BinRow_acc> load_lee_rows_acc(const std::string& lee_csv_path,
                                                 std::unordered_map<int,size_t>& bin_to_index) {
    std::ifstream fin(lee_csv_path);
    if (!fin.is_open()) {
        fatal_acc("Cannot open Lee CSV: " + lee_csv_path);
    }

    std::string header_line;
    if (!std::getline(fin, header_line)) {
        fatal_acc("Lee CSV appears empty: " + lee_csv_path);
    }
    std::vector<std::string> header = split_csv_line_acc(header_line);

    LeeCsvCols_acc cols_lee = detect_lee_columns_acc(header);
    auto H = build_header_index_acc(header);

    const std::vector<std::string> required = {
        "valid bin",
        "acceptance, ep->epg, sim(KM15)"
    };
    require_columns_acc(H, required, "Lee CSV");

    std::vector<BinRow_acc> rows;
    rows.reserve(3000);

    std::string line;
    int input_rows = 0;
    int kept_rows  = 0;

    while (std::getline(fin, line)) {
        ++input_rows;
        if (line.empty()) continue;

        std::vector<std::string> cols = split_csv_line_acc(line);

        const std::string& valid_s = get_col_ref_acc(cols, H, "valid bin");
        int valid = ToInt_acc(valid_s);
        if (valid != 1) continue;

        if ((int)cols.size() <= cols_lee.c_phi_avg) {
            fatal_acc("Row in Lee CSV has fewer columns than expected based on header detection.");
        }

        BinRow_acc r;
        r.bin_index   = ToInt_acc(cols[cols_lee.c_bin]);
        r.xBmin       = ToDouble_acc(cols[cols_lee.c_xb_min]);
        r.xBmax       = ToDouble_acc(cols[cols_lee.c_xb_max]);
        r.Q2min       = ToDouble_acc(cols[cols_lee.c_q2_min]);
        r.Q2max       = ToDouble_acc(cols[cols_lee.c_q2_max]);
        r.tmin        = ToDouble_acc(cols[cols_lee.c_t_min]);
        r.tmax        = ToDouble_acc(cols[cols_lee.c_t_max]);
        r.phiavg      = ToDouble_acc(cols[cols_lee.c_phi_avg]);

        double lee_raw = ToDouble_acc(
            get_col_ref_acc(cols, H, "acceptance, ep->epg, sim(KM15)")
        );
        r.lee_acc         = lee_raw;
        r.my_acc_inb      = 0.0;
        r.my_acc_inb_err  = 0.0;
        r.my_acc_out      = 0.0;
        r.my_acc_out_err  = 0.0;

        if (bin_to_index.find(r.bin_index) != bin_to_index.end()) {
            fatal_acc("Duplicate bin index in Lee CSV: " + std::to_string(r.bin_index));
        }

        bin_to_index[r.bin_index] = rows.size();
        rows.push_back(r);
        ++kept_rows;
    }

    info_acc("Lee CSV rows read: " + std::to_string(input_rows));
    info_acc("Lee valid rows kept (valid bin == 1): " + std::to_string(kept_rows));
    return rows;
}

static void fill_hayward_acc(const std::string& hayward_csv_path,
                             const std::unordered_map<int,size_t>& bin_to_index,
                             std::vector<BinRow_acc>& rows) {
    std::ifstream fin(hayward_csv_path);
    if (!fin.is_open()) {
        fatal_acc("Cannot open Hayward CSV: " + hayward_csv_path);
    }

    std::string header_line;
    if (!std::getline(fin, header_line)) {
        fatal_acc("Hayward CSV appears empty: " + hayward_csv_path);
    }
    std::vector<std::string> header = split_csv_line_acc(header_line);
    auto H = build_header_index_acc(header);

    const std::vector<std::string> required = {
        "bin index",
        "valid bin",
        "acceptance, Fa18 Inb",
        "acceptance, Fa18 Out"
    };
    require_columns_acc(H, required, "Hayward CSV");

    std::string line;
    int input_rows = 0;
    int matched    = 0;

    while (std::getline(fin, line)) {
        ++input_rows;
        if (line.empty()) continue;

        std::vector<std::string> cols = split_csv_line_acc(line);
        const std::string& valid_s = get_col_ref_acc(cols, H, "valid bin");
        int valid = ToInt_acc(valid_s);
        if (valid != 1) continue;

        int bin_index = ToInt_acc(get_col_ref_acc(cols, H, "bin index"));
        auto it = bin_to_index.find(bin_index);
        if (it == bin_to_index.end()) {
            continue;
        }

        AccValue_acc av_inb = parse_acc_triplet_acc(
            get_col_ref_acc(cols, H, "acceptance, Fa18 Inb")
        );
        AccValue_acc av_out = parse_acc_triplet_acc(
            get_col_ref_acc(cols, H, "acceptance, Fa18 Out")
        );

        BinRow_acc& r = rows[it->second];
        if (av_inb.valid) {
            r.my_acc_inb     = av_inb.value;
            r.my_acc_inb_err = av_inb.stat_err;
        }
        if (av_out.valid) {
            r.my_acc_out     = av_out.value;
            r.my_acc_out_err = av_out.stat_err;
        }
        ++matched;
    }

    info_acc("Hayward CSV rows read: " + std::to_string(input_rows));
    info_acc("Hayward valid rows matched to Lee bins: " + std::to_string(matched));
}

// ---------- driver ----------

void plot_acceptance_cross_checks(const std::string& lee_csv_path,
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
    auto rows = load_lee_rows_acc(lee_csv_path, bin_to_index);

    // 2) Load Hayward CSV
    fill_hayward_acc(hayward_csv_path, bin_to_index, rows);

    // 3) Build axis sets and per-panel maps
    AxisSets_acc ax = build_axes_from_rows_acc(rows);
    PerPanel_acc pp = map_to_panels_acc(rows, ax);

    info_acc("Axis xB bins: " + std::to_string(ax.xB.size()));

    auto make_fetchAll = [&](int ix) {
        return [&, ix](int iQcol, int irow,
                       PanelData_acc& inb,
                       PanelData_acc& out,
                       PanelData_acc& lee) {
            auto key = std::make_tuple(ix, iQcol, irow);

            auto itInb = pp.inb.find(key);
            if (itInb != pp.inb.end()) {
                inb = itInb->second;
            } else {
                inb = PanelData_acc();
            }

            auto itOut = pp.out.find(key);
            if (itOut != pp.out.end()) {
                out = itOut->second;
            } else {
                out = PanelData_acc();
            }

            auto itLee = pp.lee.find(key);
            if (itLee != pp.lee.end()) {
                lee = itLee->second;
            } else {
                lee = PanelData_acc();
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
            Form("Acceptance: 10.6 GeV   x_{B} #in [%.3g, %.3g]", xb_lo, xb_hi);
        const std::string title_ratio  =
            Form("Acceptance ratio (Hayward/Lee): 10.6 GeV   x_{B} #in [%.3g, %.3g]",
                 xb_lo, xb_hi);

        auto fetchAll = make_fetchAll(ix);

        const std::string f_counts =
            (fs::path(output_base_dir) / Form("acceptance_counts_xB_%d.png", ix)).string();
        const std::string f_ratio  =
            (fs::path(output_base_dir) / Form("acceptance_ratio_xB_%d.png",  ix)).string();

        draw_one_canvas_acc(title_counts, Q2s, Ts, fetchAll, f_counts, /*draw_ratio_only=*/false);
        draw_one_canvas_acc(title_ratio,  Q2s, Ts, fetchAll, f_ratio,   /*draw_ratio_only=*/true);

        info_acc("Saved: " + f_counts);
        info_acc("Saved: " + f_ratio);
    }
}