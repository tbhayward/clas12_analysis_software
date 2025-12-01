// rad_correction_cross_check.cpp
// -----------------------------------------------------------------------------
// Cross-check of radiative correction factors Frad (Hayward pass-2 vs Lee pass-1)
// using *only* CSVs.
// - Lee CSV (pass-1): e.g. imports/all_bin_v3.csv
// - Hayward CSV (pass-2): e.g. output/csvs/dvcs_pass2_analysis.csv
//
// For each valid bin, we compare:
//
//   Lee:     "Frad" (single value)
//   Hayward: "Frad, 10.6 GeV" (three-tuple: value, stat_err, syst_err)
//
// Since radiative corrections depend on beam energy (not magnet polarity),
// there is no inbending/outbending split - just one overall comparison.
//
// Output filenames:
//   frad_counts_xB_<ix>.png
//   frad_ratio_xB_<ix>.png
//
// -----------------------------------------------------------------------------

#include "rad_correction_cross_check.h"

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

static inline void info(const std::string& s) {
    std::cout << "[frad] " << s << std::endl;
}

static inline void warn(const std::string& s) {
    std::cout << "[frad][warn] " << s << std::endl;
}

static inline void fatal(const std::string& s) {
    std::cerr << "[frad][FATAL] " << s << std::endl;
    throw std::runtime_error(s);
}

static inline std::string slower(std::string s) {
    for (auto& c : s) {
        c = (char)std::tolower((unsigned char)c);
    }
    return s;
}

// trim leading/trailing whitespace
static inline std::string trim(const std::string& s) {
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

static std::vector<std::string> split_csv_line(const std::string& line) {
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
build_header_index(const std::vector<std::string>& header) {
    std::unordered_map<std::string, int> m;
    for (int i = 0; i < (int)header.size(); ++i) {
        m[header[i]] = i;
    }
    return m;
}

static const std::string& get_col_ref(const std::vector<std::string>& row,
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

static double ToDouble(const std::string& s) {
    if (s.empty()) return 0.0;
    char* endp = nullptr;
    double v = std::strtod(s.c_str(), &endp);
    if (endp == s.c_str()) return 0.0;
    return v;
}

static int ToInt(const std::string& s) {
    if (s.empty()) return 0;
    char* endp = nullptr;
    long v = std::strtol(s.c_str(), &endp, 10);
    if (endp == s.c_str()) return 0;
    return (int)v;
}

static void require_columns(const std::unordered_map<std::string,int>& idx,
                            const std::vector<std::string>& cols,
                            const std::string& which_csv) {
    std::vector<std::string> missing;
    for (const auto& c : cols) {
        if (idx.find(c) == idx.end()) {
            missing.push_back(c);
        }
    }
    if (!missing.empty()) {
        std::ostringstream oss;
        oss << "Missing columns in " << which_csv << ": ";
        for (size_t i = 0; i < missing.size(); ++i) {
            if (i > 0) oss << ", ";
            oss << '"' << missing[i] << '"';
        }
        fatal(oss.str());
    }
}

// ---------- Three-tuple parser for Hayward Frad column ----------
// Format: "(value, stat_err, syst_err)" - we extract value and stat_err

struct FradValue {
    double value    = 0.0;
    double stat_err = 0.0;
    bool   valid    = false;
};

static FradValue parse_frad_triplet(const std::string& s) {
    FradValue fv;
    if (s.empty()) return fv;

    size_t p0 = s.find('(');
    size_t p1 = s.find(')');
    if (p0 == std::string::npos || p1 == std::string::npos || p1 <= p0) {
        return fv;
    }

    std::string inner = s.substr(p0 + 1, p1 - p0 - 1);

    std::vector<std::string> parts;
    std::string cur;
    for (char c : inner) {
        if (c == ',') {
            parts.push_back(trim(cur));
            cur.clear();
        } else {
            cur.push_back(c);
        }
    }
    parts.push_back(trim(cur));

    if (parts.size() >= 2) {
        fv.value    = ToDouble(parts[0]);
        fv.stat_err = ToDouble(parts[1]);
        fv.valid    = true;
    }

    return fv;
}

// ---------- Lee CSV column detection ----------

struct LeeCsvCols {
    int c_bin;
    int c_xb_min;
    int c_xb_max;
    int c_q2_min;
    int c_q2_max;
    int c_t_min;
    int c_t_max;
    int c_phi_avg;
};

static int find_col_alias(const std::vector<std::string>& header,
                          const std::vector<std::string>& names) {
    for (size_t i = 0; i < header.size(); ++i) {
        std::string h = slower(trim(header[i]));
        for (const auto& raw_name : names) {
            std::string n = slower(trim(raw_name));
            if (!n.empty() && h == n) {
                return (int)i;
            }
        }
    }
    return -1;
}

static LeeCsvCols detect_lee_columns(const std::vector<std::string>& header) {
    if (header.empty()) {
        fatal("Empty header row in Lee CSV");
    }

    LeeCsvCols cols;
    cols.c_bin     = -1;
    cols.c_xb_min  = -1;
    cols.c_xb_max  = -1;
    cols.c_q2_min  = -1;
    cols.c_q2_max  = -1;
    cols.c_t_min   = -1;
    cols.c_t_max   = -1;
    cols.c_phi_avg = -1;

    // 1) Bin index: try named column first, then unlabeled first column
    cols.c_bin = find_col_alias(header, { "bin index", "bin", "idx" });
    if (cols.c_bin < 0) {
        std::string h0 = trim(header[0]);
        if (h0.empty()) {
            cols.c_bin = 0;
        } else {
            std::ostringstream oss;
            oss << "Could not locate bin index column in Lee CSV. "
                << "Tried names: \"bin index\", \"bin\", \"idx\" and unlabeled first column.\n"
                << "Header[0] is \"" << header[0] << "\".";
            fatal(oss.str());
        }
    }

    // 2) xB, Q2, phi
    cols.c_xb_min  = find_col_alias(header, { "xBmin", "xbmin", "xB_min", "xb_min" });
    cols.c_xb_max  = find_col_alias(header, { "xBmax", "xbmax", "xB_max", "xb_max" });
    cols.c_q2_min  = find_col_alias(header, { "Q2min", "q2min", "Q2_min", "q2_min" });
    cols.c_q2_max  = find_col_alias(header, { "Q2max", "q2max", "Q2_max", "q2_max" });
    cols.c_phi_avg = find_col_alias(header, { "phiavg", "phi_avg", "phi_average" });

    // 3) |t| min and max
    cols.c_t_min = find_col_alias(header, { "t_abs_min", "tmin", "t_min" });
    cols.c_t_max = find_col_alias(header, { "t_abs_max", "tmax", "t_max" });

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
        fatal(oss.str());
    }

    return cols;
}

// ---------- bin / axis structs ----------

struct AxisSets {
    std::vector<std::pair<double,double>> xB;
    std::map<int, std::vector<std::pair<double,double>>> Q2_by_ix;
    std::map<int, std::vector<std::pair<double,double>>> t_by_ix;
};

struct BinRow {
    int    bin_index   = 0;
    double xBmin       = 0.0;
    double xBmax       = 0.0;
    double Q2min       = 0.0;
    double Q2max       = 0.0;
    double tmin        = 0.0;
    double tmax        = 0.0;
    double phiavg      = 0.0;

    double lee_frad    = 0.0;  // Lee Frad (single value, no error)
    double my_frad     = 0.0;  // Hayward Frad, 10.6 GeV (value)
    double my_frad_err = 0.0;  // Hayward Frad stat error
};

static AxisSets build_axes_from_rows(const std::vector<BinRow>& rows) {
    std::set<std::pair<double,double>> xbset;
    std::map<std::pair<double,double>, std::set<std::pair<double,double>>> q2set_by_xb;
    std::map<std::pair<double,double>, std::set<std::pair<double,double>>> tset_by_xb;

    for (const auto& r : rows) {
        const auto xb = std::make_pair(r.xBmin, r.xBmax);
        xbset.insert(xb);
        q2set_by_xb[xb].insert({r.Q2min, r.Q2max});
        tset_by_xb[xb].insert({r.tmin,  r.tmax});
    }

    AxisSets ax;
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

static inline int find_index(const std::pair<double,double>& r,
                             const std::vector<std::pair<double,double>>& v) {
    for (int i = 0; i < (int)v.size(); ++i) {
        if (v[i] == r) return i;
    }
    return -1;
}

// Per-panel data
struct PanelData {
    std::vector<double> phi;
    std::vector<double> val;
    std::vector<double> err;
};

struct PerPanel {
    // key = (ix, iQ, it)
    std::map<std::tuple<int,int,int>, PanelData> lee;
    std::map<std::tuple<int,int,int>, PanelData> hayward;
};

static PerPanel map_to_panels(const std::vector<BinRow>& rows,
                              const AxisSets& ax) {
    PerPanel pp;

    for (const auto& r : rows) {
        const auto xb = std::make_pair(r.xBmin, r.xBmax);
        const int ix  = find_index(xb, ax.xB);
        if (ix < 0) continue;

        const auto& Q2s = ax.Q2_by_ix.at(ix);
        const auto& Ts  = ax.t_by_ix.at(ix);

        const int iQ = find_index({r.Q2min, r.Q2max}, Q2s);
        const int it = find_index({r.tmin,   r.tmax},  Ts);
        if (iQ < 0 || it < 0) continue;

        auto key = std::make_tuple(ix, iQ, it);

        // Lee data (no error bars)
        if (r.lee_frad > 0.0) {
            pp.lee[key].phi.push_back(r.phiavg);
            pp.lee[key].val.push_back(r.lee_frad);
            pp.lee[key].err.push_back(0.0);
        }

        // Hayward data (with error bars)
        if (r.my_frad > 0.0) {
            pp.hayward[key].phi.push_back(r.phiavg);
            pp.hayward[key].val.push_back(r.my_frad);
            pp.hayward[key].err.push_back(r.my_frad_err);
        }
    }

    // Sort each panel's data by phi
    auto sort_panel = [](PanelData& pd) {
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

static inline void degreeTicks(double xmin, double ymin, double xmax, double labelSize) {
    TGaxis* ax = new TGaxis(xmin, ymin, xmax, ymin, 0.0, 360.0, 4, "");
    ax->SetLabelFont(42);
    ax->SetLabelSize(labelSize);
    ax->SetLabelOffset(0.012);
    ax->SetTitle("");
    ax->SetTickSize(0.02);
    ax->Draw();
}

static void draw_degree_ticks_here(double labelSize) {
    degreeTicks(gPad->GetUxmin(), gPad->GetUymin(), gPad->GetUxmax(), labelSize);
}

// Markers + error bars
static TGraphErrors* graph_pe1(const std::vector<double>& X,
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

static std::string safe_canvas_name(const std::string& out_png) {
    return fs::path(out_png).filename().string();
}

static void draw_one_canvas(const std::string& title,
                            const std::vector<std::pair<double,double>>& Q2s,
                            const std::vector<std::pair<double,double>>& Ts,
                            const std::function<void(int,int,PanelData&,PanelData&)>& fetchBoth,
                            const std::string& out_png,
                            bool draw_ratio_only) {
    const int nrows = (int)Ts.size();
    const int ncols = (int)Q2s.size();
    if (nrows == 0 || ncols == 0) return;

    const int W = 320 * ncols + 220;
    const int H = 260 * nrows + 260;

    const std::string cname = safe_canvas_name(out_png);
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

            PanelData hayward, lee;
            fetchBoth(ccol, r, hayward, lee);

            // Fixed y-axis range
            // Frad is typically close to 1, so use 0.8 to 1.2 for values, 0.8 to 1.2 for ratio
            double ymin = draw_ratio_only ? 0.8 : 0.8;
            double ymax = draw_ratio_only ? 1.2 : 1.2;

            TH1* frame = gPad->DrawFrame(0.0, ymin, 360.0, ymax);
            frame->GetXaxis()->SetTitle("#phi (deg)");
            frame->GetYaxis()->SetTitle(draw_ratio_only ? "Hayward / Lee" : "F_{rad}");
            frame->GetXaxis()->CenterTitle();
            frame->GetYaxis()->CenterTitle();
            frame->GetXaxis()->SetNdivisions(505);
            frame->GetXaxis()->SetTitleSize(0.060);
            frame->GetYaxis()->SetTitleSize(0.060);
            frame->GetYaxis()->SetLabelSize(0.050);
            frame->GetXaxis()->SetLabelSize(0.050);
            frame->GetXaxis()->SetTitleOffset(1.10);
            frame->GetYaxis()->SetTitleOffset(1.55);

            draw_degree_ticks_here(0.050);

            TLatex lab;
            lab.SetNDC();
            lab.SetTextSize(0.070);
            lab.SetTextAlign(11);
            lab.SetTextFont(42);
            lab.DrawLatex(0.15, 0.92,
                Form("Q^{2} #in [%.2g, %.2g], -t #in [%.2g, %.2g]",
                     Q2s[ccol].first, Q2s[ccol].second,
                     Ts[r].first,     Ts[r].second));

            if (draw_ratio_only) {
                // Compute ratios by matching phi positions
                const double tol = 20.0;
                std::vector<double> x, y, ey;

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
                        double R = (hayward.val[i] <= 0.0) ? 0.0 : hayward.val[i] / lee.val[jbest];
                        double eR = 0.0;
                        if (hayward.val[i] > 0.0 && hayward.err[i] > 0.0) {
                            eR = R * (hayward.err[i] / hayward.val[i]);
                        }
                        x.push_back(hayward.phi[i]);
                        y.push_back(R);
                        ey.push_back(eR);
                    }
                }

                graph_pe1(x, y, ey, 20, black);

                TLine* one = new TLine(0.0, 1.0, 360.0, 1.0);
                one->SetLineStyle(2);
                one->SetLineWidth(2);
                one->SetLineColor(orange);
                one->Draw("SAME");
            } else {
                // Plot both series
                graph_pe1(hayward.phi, hayward.val, hayward.err, 20, black);  // Hayward with errors
                graph_pe1(lee.phi, lee.val, lee.err, 24, orange);             // Lee (errors are 0)
            }
        }
    }

    c->SaveAs(out_png.c_str());

    delete legTop;
    delete c;
}

// ---------- CSV loaders ----------

static std::vector<BinRow> load_lee_rows(const std::string& lee_csv_path,
                                         std::unordered_map<int,size_t>& bin_to_index) {
    std::ifstream fin(lee_csv_path);
    if (!fin.is_open()) {
        fatal("Cannot open Lee CSV: " + lee_csv_path);
    }

    std::string header_line;
    if (!std::getline(fin, header_line)) {
        fatal("Lee CSV appears empty: " + lee_csv_path);
    }
    std::vector<std::string> header = split_csv_line(header_line);

    LeeCsvCols cols_lee = detect_lee_columns(header);
    auto H = build_header_index(header);

    const std::vector<std::string> required = {
        "valid bin",
        "Frad"
    };
    require_columns(H, required, "Lee CSV");

    std::vector<BinRow> rows;
    rows.reserve(3000);

    std::string line;
    int input_rows = 0;
    int kept_rows  = 0;

    while (std::getline(fin, line)) {
        ++input_rows;
        if (line.empty()) continue;

        std::vector<std::string> cols = split_csv_line(line);

        const std::string& valid_s = get_col_ref(cols, H, "valid bin");
        int valid = ToInt(valid_s);
        if (valid != 1) continue;

        if ((int)cols.size() <= cols_lee.c_phi_avg) {
            fatal("Row in Lee CSV has fewer columns than expected based on header detection.");
        }

        BinRow r;
        r.bin_index   = ToInt(cols[cols_lee.c_bin]);
        r.xBmin       = ToDouble(cols[cols_lee.c_xb_min]);
        r.xBmax       = ToDouble(cols[cols_lee.c_xb_max]);
        r.Q2min       = ToDouble(cols[cols_lee.c_q2_min]);
        r.Q2max       = ToDouble(cols[cols_lee.c_q2_max]);
        r.tmin        = ToDouble(cols[cols_lee.c_t_min]);
        r.tmax        = ToDouble(cols[cols_lee.c_t_max]);
        r.phiavg      = ToDouble(cols[cols_lee.c_phi_avg]);

        r.lee_frad    = ToDouble(get_col_ref(cols, H, "Frad"));
        r.my_frad     = 0.0;
        r.my_frad_err = 0.0;

        if (bin_to_index.find(r.bin_index) != bin_to_index.end()) {
            fatal("Duplicate bin index in Lee CSV: " + std::to_string(r.bin_index));
        }

        bin_to_index[r.bin_index] = rows.size();
        rows.push_back(r);
        ++kept_rows;
    }

    info("Lee CSV rows read: " + std::to_string(input_rows));
    info("Lee valid rows kept (valid bin == 1): " + std::to_string(kept_rows));
    return rows;
}

static void fill_hayward_frad(const std::string& hayward_csv_path,
                              const std::unordered_map<int,size_t>& bin_to_index,
                              std::vector<BinRow>& rows) {
    std::ifstream fin(hayward_csv_path);
    if (!fin.is_open()) {
        fatal("Cannot open Hayward CSV: " + hayward_csv_path);
    }

    std::string header_line;
    if (!std::getline(fin, header_line)) {
        fatal("Hayward CSV appears empty: " + hayward_csv_path);
    }
    std::vector<std::string> header = split_csv_line(header_line);
    auto H = build_header_index(header);

    const std::vector<std::string> required = {
        "bin index",
        "valid bin",
        "Frad, 10.6 GeV"
    };
    require_columns(H, required, "Hayward CSV");

    std::string line;
    int input_rows = 0;
    int matched    = 0;

    while (std::getline(fin, line)) {
        ++input_rows;
        if (line.empty()) continue;

        std::vector<std::string> cols = split_csv_line(line);
        const std::string& valid_s = get_col_ref(cols, H, "valid bin");
        int valid = ToInt(valid_s);
        if (valid != 1) continue;

        int bin_index = ToInt(get_col_ref(cols, H, "bin index"));
        auto it = bin_to_index.find(bin_index);
        if (it == bin_to_index.end()) {
            continue;
        }

        FradValue fv = parse_frad_triplet(get_col_ref(cols, H, "Frad, 10.6 GeV"));

        BinRow& r = rows[it->second];
        if (fv.valid) {
            r.my_frad     = fv.value;
            r.my_frad_err = fv.stat_err;
        }
        ++matched;
    }

    info("Hayward CSV rows read: " + std::to_string(input_rows));
    info("Hayward valid rows matched to Lee bins: " + std::to_string(matched));
}

// ---------- driver ----------

void plot_rad_correction_cross_checks(const std::string& lee_csv_path,
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
    auto rows = load_lee_rows(lee_csv_path, bin_to_index);

    // 2) Load Hayward CSV
    fill_hayward_frad(hayward_csv_path, bin_to_index, rows);

    // 3) Build axis sets and per-panel maps
    AxisSets ax = build_axes_from_rows(rows);
    PerPanel pp = map_to_panels(rows, ax);

    info("Axis xB bins: " + std::to_string(ax.xB.size()));

    auto make_fetchBoth = [&](int ix) {
        return [&, ix](int iQcol, int irow,
                       PanelData& hayward,
                       PanelData& lee) {
            auto key = std::make_tuple(ix, iQcol, irow);

            auto itH = pp.hayward.find(key);
            if (itH != pp.hayward.end()) {
                hayward = itH->second;
            } else {
                hayward = PanelData();
            }

            auto itL = pp.lee.find(key);
            if (itL != pp.lee.end()) {
                lee = itL->second;
            } else {
                lee = PanelData();
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
            Form("Radiative correction F_{rad}: 10.6 GeV   x_{B} #in [%.3g, %.3g]", xb_lo, xb_hi);
        const std::string title_ratio  =
            Form("F_{rad} ratio (Hayward/Lee): 10.6 GeV   x_{B} #in [%.3g, %.3g]", xb_lo, xb_hi);

        auto fetchBoth = make_fetchBoth(ix);

        const std::string f_counts =
            (fs::path(output_base_dir) / Form("frad_counts_xB_%d.png", ix)).string();
        const std::string f_ratio  =
            (fs::path(output_base_dir) / Form("frad_ratio_xB_%d.png",  ix)).string();

        draw_one_canvas(title_counts, Q2s, Ts, fetchBoth, f_counts, /*draw_ratio_only=*/false);
        draw_one_canvas(title_ratio,  Q2s, Ts, fetchBoth, f_ratio,   /*draw_ratio_only=*/true);

        info("Saved: " + f_counts);
        info("Saved: " + f_ratio);
    }
}