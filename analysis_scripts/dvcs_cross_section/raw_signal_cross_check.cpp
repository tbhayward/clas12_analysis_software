// raw_signal_cross_check.cpp
// -----------------------------------------------------------------------------
// Cross-check of raw yields (Hayward pass-2 vs Lee pass-1) using *only* CSVs.
// - Lee CSV (pass-1): e.g. imports/all_bin_v3.csv
// - Hayward CSV (pass-2): e.g. output/csvs/dvcs_pass2_analysis.csv
//
// For each valid bin, we:
//
//   Lee Fa18 Inb  = sum of
//      "raw yield, ep->epg, (FD, FD), exp, inbending"
//      "raw yield, ep->epg, (CD, FD), exp, inbending"
//      "raw yield, ep->epg, (CD, FT), exp, inbending"
//
//   Lee Fa18 Out  = sum of
//      "raw yield, ep->epg, (FD, FD), exp, outbending"
//      "raw yield, ep->epg, (CD, FD), exp, outbending"
//      "raw yield, ep->epg, (CD, FT), exp, outbending"
//
//   Hayward Fa18 Inb (unpol) = sum of
//      "raw yield, ep->epg, (FD, FD), exp, Fa18 Inb, unpol"
//      "raw yield, ep->epg, (CD, FD), exp, Fa18 Inb, unpol"
//      "raw yield, ep->epg, (CD, FT), exp, Fa18 Inb, unpol"
//
//   Hayward Fa18 Out (unpol) = sum of
//      "raw yield, ep->epg, (FD, FD), exp, Fa18 Out, unpol"
//      "raw yield, ep->epg, (CD, FD), exp, Fa18 Out, unpol"
//      "raw yield, ep->epg, (CD, FT), exp, Fa18 Out, unpol"
//
// We then build axis sets in (xB, Q^2, -t). For each (xB, Q^2, -t) cell,
// we take *all* rows matching those ranges and use their provided phiavg
// values directly as the x-coordinates. We do NOT rebin phi into 12
// uniform bins.
//
// Finally we produce, for each xB bin and for Fa18 Inb / Fa18 Out:
//
//   1) Raw counts vs phiavg: Hayward (black) vs Lee (orange)
//   2) Ratio Hayward/Lee vs phiavg with Poisson errors and a y=1 line
//
// Output filenames:
//   raw_counts_fa18_inb_xB_<ix>.png
//   raw_ratio_fa18_inb_xB_<ix>.png
//   raw_counts_fa18_out_xB_<ix>.png
//   raw_ratio_fa18_out_xB_<ix>.png
//
// -----------------------------------------------------------------------------

#include "raw_signal_cross_check.h"

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
    std::cout << "[cross] " << s << std::endl;
}

static inline void warn(const std::string& s) {
    std::cout << "[cross][warn] " << s << std::endl;
}

static inline void fatal(const std::string& s) {
    std::cerr << "[cross][FATAL] " << s << std::endl;
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

// ---------- Lee CSV column aliases and detection ----------

struct LeeCsvCols {
    int c_bin;
    int c_xb_min;
    int c_xb_max;
    int c_q2_min;
    int c_q2_max;
    int c_t_min;
    int c_t_max;
    int c_phi_min;
    int c_phi_max;
    int c_phiavg;
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
        // endfor
    }
    // endfor
    return -1;
}

static LeeCsvCols detect_lee_columns(const std::vector<std::string>& header) {
    if (header.empty()) {
        fatal("[cross][FATAL] Empty header row in Lee CSV");
    }

    LeeCsvCols cols;
    cols.c_bin     = -1;
    cols.c_xb_min  = -1;
    cols.c_xb_max  = -1;
    cols.c_q2_min  = -1;
    cols.c_q2_max  = -1;
    cols.c_t_min   = -1;
    cols.c_t_max   = -1;
    cols.c_phi_min = -1;
    cols.c_phi_max = -1;
    cols.c_phiavg  = -1;

    // 1) Bin index: try named column first, then unlabeled first column
    cols.c_bin = find_col_alias(header, { "bin index", "bin", "idx" });
    if (cols.c_bin < 0) {
        std::string h0 = trim(header[0]);
        if (h0.empty()) {
            // Lee pass-1: first column is bin index, header cell is blank
            cols.c_bin = 0;
        } else {
            std::ostringstream oss;
            oss << "[cross][FATAL] Could not locate bin index column in Lee CSV. "
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
    cols.c_phi_min = find_col_alias(header, { "phimin", "phi_min", "phi_minimum" });
    cols.c_phi_max = find_col_alias(header, { "phimax", "phi_max", "phi_maximum" });
    cols.c_phiavg  = find_col_alias(header, { "phiavg", "phi_avg", "phi_average" });

    // 3) |t| min and max: accept both our names and Lee's pass-1 names.
    cols.c_t_min   = find_col_alias(header, { "tmin", "t_min", "t_abs_min" });
    cols.c_t_max   = find_col_alias(header, { "tmax", "t_max", "t_abs_max" });

    // 4) Validate required columns and give a clear fatal if something is missing.
    std::vector<std::string> missing;

    if (cols.c_bin     < 0) missing.push_back("bin index");
    if (cols.c_xb_min  < 0) missing.push_back("xBmin");
    if (cols.c_xb_max  < 0) missing.push_back("xBmax");
    if (cols.c_q2_min  < 0) missing.push_back("Q2min");
    if (cols.c_q2_max  < 0) missing.push_back("Q2max");
    if (cols.c_t_min   < 0) missing.push_back("tmin/t_abs_min");
    if (cols.c_t_max   < 0) missing.push_back("tmax/t_abs_max");
    if (cols.c_phi_min < 0) missing.push_back("phimin");
    if (cols.c_phi_max < 0) missing.push_back("phimax");
    if (cols.c_phiavg  < 0) missing.push_back("phiavg");

    if (!missing.empty()) {
        std::ostringstream oss;
        oss << "[cross][FATAL] Missing required columns in Lee CSV: ";
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

// Row struct combining bin definition and both analyses' counts
struct BinRow {
    int    bin_index = 0;
    double xBmin     = 0.0;
    double xBmax     = 0.0;
    double Q2min     = 0.0;
    double Q2max     = 0.0;
    double tmin      = 0.0;
    double tmax      = 0.0;
    double phiavg    = 0.0;

    double lee_inb   = 0.0; // Lee Fa18 Inb (inbending)
    double lee_out   = 0.0; // Lee Fa18 Out (outbending)
    double my_inb    = 0.0; // Hayward Fa18 Inb, unpol
    double my_out    = 0.0; // Hayward Fa18 Out, unpol
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

// Markers + error bars only (no connecting line)
static TGraphErrors* graph_pe1(const std::vector<double>& X,
                               const std::vector<double>& Y,
                               const std::vector<double>& EY,
                               int markerStyle, int color) {
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

// Fill function type: for given (Q2 index, t index) returns variable-length
// vectors of phiavg, Hayward counts, and Lee counts.
using FillFunc = std::function<void(int,int,
                                    std::vector<double>&,
                                    std::vector<double>&,
                                    std::vector<double>&)>;

static double compute_canvas_ymax(
    bool ratio_mode,
    const std::vector<std::pair<double,double>>& Q2s,
    const std::vector<std::pair<double,double>>& Ts,
    const FillFunc& fillBoth
) {
    double ymax = 0.0;
    const int nrows = (int)Ts.size();
    const int ncols = (int)Q2s.size();

    for (int r = 0; r < nrows; ++r) {
        for (int ccol = 0; ccol < ncols; ++ccol) {
            std::vector<double> phi, A, B;
            fillBoth(ccol, r, phi, A, B);

            if (!ratio_mode) {
                for (size_t i = 0; i < A.size(); ++i) {
                    ymax = std::max(ymax, A[i]);
                }
                for (size_t i = 0; i < B.size(); ++i) {
                    ymax = std::max(ymax, B[i]);
                }
            } else {
                const size_t n = std::min(A.size(), B.size());
                for (size_t i = 0; i < n; ++i) {
                    const double NL = B[i];
                    const double NH = A[i];
                    if (NL <= 0.0) continue;
                    const double R = (NH <= 0.0) ? 0.0 : NH / NL;
                    double eR = 0.0;
                    if (NH > 0.0) {
                        eR = R * std::sqrt(1.0 / NH + 1.0 / NL);
                    }
                    ymax = std::max(ymax, R + eR);
                }
            }
        }
    }

    if (ymax <= 0.0) ymax = 1.0;
    if (ratio_mode)  ymax = std::max(ymax, 1.0);
    return ymax;
}

static void draw_one_canvas(const std::string& title,
                            const std::vector<std::pair<double,double>>& Q2s,
                            const std::vector<std::pair<double,double>>& Ts,
                            const FillFunc& fillBoth,
                            const std::string& out_png,
                            bool draw_ratio_only) {
    const int nrows = (int)Ts.size();
    const int ncols = (int)Q2s.size();
    if (nrows == 0 || ncols == 0) return;

    const double canvas_ymax = compute_canvas_ymax(draw_ratio_only, Q2s, Ts, fillBoth);

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

            std::vector<double> phi, A, B;
            fillBoth(ccol, r, phi, A, B);

            TH1* frame = gPad->DrawFrame(0.0, 0.0, 360.0, canvas_ymax);
            frame->GetXaxis()->SetTitle("#phi (deg)");
            frame->GetYaxis()->SetTitle(draw_ratio_only ? "Hayward / Lee" : "Raw counts");
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

            if (phi.empty()) {
                // No bins for this (Q2, t) cell; leave the frame empty.
                continue;
            }

            if (draw_ratio_only) {
                std::vector<double> x, y, ey;
                const size_t n = std::min(A.size(), B.size());
                x.reserve(n);
                y.reserve(n);
                ey.reserve(n);

                for (size_t i = 0; i < n; ++i) {
                    const double NL = B[i];
                    const double NH = A[i];
                    if (NL <= 0.0) continue;
                    const double R = (NH <= 0.0) ? 0.0 : NH / NL;
                    double eR = 0.0;
                    if (NH > 0.0) {
                        eR = R * std::sqrt(1.0 / NH + 1.0 / NL);
                    }
                    x.push_back(phi[i]);
                    y.push_back(R);
                    ey.push_back(eR);
                }

                if (!x.empty()) {
                    graph_pe1(x, y, ey, 20, black);

                    TLine* one = new TLine(0.0, 1.0, 360.0, 1.0);
                    one->SetLineStyle(2);
                    one->SetLineWidth(2);
                    one->SetLineColor(orange);
                    one->Draw("SAME");
                }
            } else {
                const size_t n = phi.size();
                std::vector<double> eyA(n, 0.0), eyB(n, 0.0);
                for (size_t i = 0; i < n; ++i) {
                    eyA[i] = (A[i] > 0.0) ? std::sqrt(A[i]) : 0.0;
                    eyB[i] = (B[i] > 0.0) ? std::sqrt(B[i]) : 0.0;
                }
                graph_pe1(phi, A, eyA, 20, black);   // Hayward (pass-2)
                graph_pe1(phi, B, eyB, 24, orange);  // Lee (pass-1)
            }
        }
    }

    c->SaveAs(out_png.c_str());

    delete legTop;
    delete c;
}

// ---------- CSV loaders for Lee and Hayward ----------

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

    // Detect bin / kinematic columns using aliases and Lee conventions.
    LeeCsvCols cols_lee = detect_lee_columns(header);

    // Build header index for named columns like "valid bin" and yield columns.
    auto H = build_header_index(header);

    // Only "valid bin" is required by name here; bin index and kinematic
    // ranges are handled by detect_lee_columns above.
    const std::vector<std::string> required = {
        "valid bin"
    };
    require_columns(H, required, "Lee CSV");

    const std::vector<std::string> lee_inb_cols = {
        "raw yield, ep->epg, (FD, FD), exp, inbending",
        "raw yield, ep->epg, (CD, FD), exp, inbending",
        "raw yield, ep->epg, (CD, FT), exp, inbending"
    };

    const std::vector<std::string> lee_out_cols = {
        "raw yield, ep->epg, (FD, FD), exp, outbending",
        "raw yield, ep->epg, (CD, FD), exp, outbending",
        "raw yield, ep->epg, (CD, FT), exp, outbending"
    };

    require_columns(H, lee_inb_cols, "Lee CSV");
    require_columns(H, lee_out_cols, "Lee CSV");

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

        if ((int)cols.size() <= cols_lee.c_phiavg) {
            fatal("Row in Lee CSV has fewer columns than expected based on header detection.");
        }

        BinRow r;
        r.bin_index = ToInt(cols[cols_lee.c_bin]);
        r.xBmin     = ToDouble(cols[cols_lee.c_xb_min]);
        r.xBmax     = ToDouble(cols[cols_lee.c_xb_max]);
        r.Q2min     = ToDouble(cols[cols_lee.c_q2_min]);
        r.Q2max     = ToDouble(cols[cols_lee.c_q2_max]);
        r.tmin      = ToDouble(cols[cols_lee.c_t_min]);
        r.tmax      = ToDouble(cols[cols_lee.c_t_max]);
        r.phiavg    = ToDouble(cols[cols_lee.c_phiavg]);

        double lee_inb_sum = 0.0;
        for (const auto& c : lee_inb_cols) {
            lee_inb_sum += ToDouble(get_col_ref(cols, H, c));
        }

        double lee_out_sum = 0.0;
        for (const auto& c : lee_out_cols) {
            lee_out_sum += ToDouble(get_col_ref(cols, H, c));
        }

        r.lee_inb = lee_inb_sum;
        r.lee_out = lee_out_sum;
        r.my_inb  = 0.0;
        r.my_out  = 0.0;

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

static void fill_hayward_counts(const std::string& hayward_csv_path,
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
        "bin index", "valid bin"
    };
    require_columns(H, required, "Hayward CSV");

    const std::vector<std::string> my_inb_cols = {
        "raw yield, ep->epg, (FD, FD), exp, Fa18 Inb, unpol",
        "raw yield, ep->epg, (CD, FD), exp, Fa18 Inb, unpol",
        "raw yield, ep->epg, (CD, FT), exp, Fa18 Inb, unpol"
    };

    const std::vector<std::string> my_out_cols = {
        "raw yield, ep->epg, (FD, FD), exp, Fa18 Out, unpol",
        "raw yield, ep->epg, (CD, FD), exp, Fa18 Out, unpol",
        "raw yield, ep->epg, (CD, FT), exp, Fa18 Out, unpol"
    };

    require_columns(H, my_inb_cols, "Hayward CSV");
    require_columns(H, my_out_cols, "Hayward CSV");

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
            fatal("Hayward CSV has bin index not present in Lee CSV: " + std::to_string(bin_index));
        }

        double my_inb_sum = 0.0;
        for (const auto& c : my_inb_cols) {
            my_inb_sum += ToDouble(get_col_ref(cols, H, c));
        }

        double my_out_sum = 0.0;
        for (const auto& c : my_out_cols) {
            my_out_sum += ToDouble(get_col_ref(cols, H, c));
        }

        BinRow& r = rows[it->second];
        r.my_inb  = my_inb_sum;
        r.my_out  = my_out_sum;
        ++matched;
    }

    info("Hayward CSV rows read: " + std::to_string(input_rows));
    info("Hayward valid rows matched to Lee bins: " + std::to_string(matched));
}

// ---------- helpers to gather per-(Q2,t) data ----------

static void gather_phi_rows(const std::vector<BinRow>& rows,
                            const AxisSets& ax,
                            int ix,
                            int iQ,
                            int it,
                            bool use_inb,
                            std::vector<double>& phi,
                            std::vector<double>& my_counts,
                            std::vector<double>& lee_counts) {
    phi.clear();
    my_counts.clear();
    lee_counts.clear();

    const auto xb  = ax.xB[ix];
    const auto& Q2s = ax.Q2_by_ix.at(ix);
    const auto& Ts  = ax.t_by_ix.at(ix);

    if (iQ < 0 || iQ >= (int)Q2s.size()) return;
    if (it < 0 || it >= (int)Ts.size())  return;

    const auto q2r = Q2s[iQ];
    const auto tr  = Ts[it];

    for (const auto& r : rows) {
        if (r.xBmin == xb.first && r.xBmax == xb.second &&
            r.Q2min == q2r.first && r.Q2max == q2r.second &&
            r.tmin  == tr.first  && r.tmax  == tr.second) {

            phi.push_back(r.phiavg);
            if (use_inb) {
                my_counts.push_back(r.my_inb);
                lee_counts.push_back(r.lee_inb);
            } else {
                my_counts.push_back(r.my_out);
                lee_counts.push_back(r.lee_out);
            }
        }
    }

    // Sort by phi so the graphs look nice.
    std::vector<size_t> order(phi.size());
    for (size_t i = 0; i < order.size(); ++i) {
        order[i] = i;
    }
    std::sort(order.begin(), order.end(),
              [&](size_t a, size_t b){ return phi[a] < phi[b]; });

    std::vector<double> phi_s, my_s, lee_s;
    phi_s.reserve(order.size());
    my_s.reserve(order.size());
    lee_s.reserve(order.size());

    for (size_t k = 0; k < order.size(); ++k) {
        size_t idx = order[k];
        phi_s.push_back(phi[idx]);
        my_s.push_back(my_counts[idx]);
        lee_s.push_back(lee_counts[idx]);
    }

    phi.swap(phi_s);
    my_counts.swap(my_s);
    lee_counts.swap(lee_s);
}

// ---------- driver ----------

void plot_raw_yield_cross_checks(const std::string& lee_csv_path,
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

    // 1) Load Lee CSV and build BinRow list + bin-index map
    std::unordered_map<int,size_t> bin_to_index;
    auto rows = load_lee_rows(lee_csv_path, bin_to_index);

    // 2) Load Hayward CSV and fill my_inb/my_out for matching bins
    fill_hayward_counts(hayward_csv_path, bin_to_index, rows);

    // 3) Build axis sets
    AxisSets ax = build_axes_from_rows(rows);

    info("Axis xB bins: " + std::to_string(ax.xB.size()));

    for (int ix = 0; ix < (int)ax.xB.size(); ++ix) {
        const auto& Q2s = ax.Q2_by_ix[ix];
        const auto& Ts  = ax.t_by_ix[ix];
        if (Q2s.empty() || Ts.empty()) continue;

        const double xb_lo = ax.xB[ix].first;
        const double xb_hi = ax.xB[ix].second;

        // Fa18 Inb
        {
            const std::string title_counts =
                Form("Raw yields: Fa18 Inb   x_{B} #in [%.3g, %.3g]", xb_lo, xb_hi);
            const std::string title_ratio  =
                Form("Raw yields ratio (Hayward/Lee): Fa18 Inb   x_{B} #in [%.3g, %.3g]", xb_lo, xb_hi);

            FillFunc fill_inb = [&](int iQcol, int irow,
                                    std::vector<double>& phi,
                                    std::vector<double>& A,
                                    std::vector<double>& B) {
                gather_phi_rows(rows, ax, ix, iQcol, irow, true, phi, A, B);
            };

            const std::string f_counts =
                (fs::path(output_base_dir) / Form("raw_counts_fa18_inb_xB_%d.png", ix)).string();
            const std::string f_ratio  =
                (fs::path(output_base_dir) / Form("raw_ratio_fa18_inb_xB_%d.png",  ix)).string();

            draw_one_canvas(title_counts, Q2s, Ts, fill_inb, f_counts, /*draw_ratio_only=*/false);
            draw_one_canvas(title_ratio,  Q2s, Ts, fill_inb, f_ratio,   /*draw_ratio_only=*/true);

            info("Saved: " + f_counts);
            info("Saved: " + f_ratio);
        }

        // Fa18 Out
        {
            const std::string title_counts =
                Form("Raw yields: Fa18 Out   x_{B} #in [%.3g, %.3g]", xb_lo, xb_hi);
            const std::string title_ratio  =
                Form("Raw yields ratio (Hayward/Lee): Fa18 Out   x_{B} #in [%.3g, %.3g]", xb_lo, xb_hi);

            FillFunc fill_out = [&](int iQcol, int irow,
                                    std::vector<double>& phi,
                                    std::vector<double>& A,
                                    std::vector<double>& B) {
                gather_phi_rows(rows, ax, ix, iQcol, irow, false, phi, A, B);
            };

            const std::string f_counts =
                (fs::path(output_base_dir) / Form("raw_counts_fa18_out_xB_%d.png", ix)).string();
            const std::string f_ratio  =
                (fs::path(output_base_dir) / Form("raw_ratio_fa18_out_xB_%d.png",  ix)).string();

            draw_one_canvas(title_counts, Q2s, Ts, fill_out, f_counts, /*draw_ratio_only=*/false);
            draw_one_canvas(title_ratio,  Q2s, Ts, fill_out, f_ratio,   /*draw_ratio_only=*/true);

            info("Saved: " + f_counts);
            info("Saved: " + f_ratio);
        }
    }
}