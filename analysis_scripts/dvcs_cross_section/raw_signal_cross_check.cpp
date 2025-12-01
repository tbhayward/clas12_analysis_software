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
// We then map each row into index space (ix, iQ, it, ip), with:
//   ix  : xB bin index (unique [xBmin,xBmax])
//   iQ  : Q^2 bin index within ix (unique [Q2min,Q2max])
//   it  : -t bin index within ix (unique [tmin,tmax])
//   ip  : phi bin index from phiavg into 12 uniform bins
//
// Finally we produce, for each xB bin and for Fa18 Inb / Fa18 Out:
//
//   1) Raw counts vs phi: Hayward (black) vs Lee (orange)
//   2) Ratio Hayward/Lee vs phi with Poisson errors and a y=1 line
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

// ---------- bin / axis structs ----------

static constexpr int N_PHI_BINS = 12;

static inline int phiBinFromDeg(double phi_deg) {
    if (!std::isfinite(phi_deg)) return -1;
    double p = std::fmod(phi_deg, 360.0);
    if (p < 0.0) p += 360.0;
    const double w = 360.0 / double(N_PHI_BINS);
    int ip = (int)std::floor(p / w);
    if (ip < 0) ip = 0;
    if (ip >= N_PHI_BINS) ip = N_PHI_BINS - 1;
    return ip;
}

static inline std::vector<double> phiCentersDeg() {
    std::vector<double> v(N_PHI_BINS);
    const double step = 360.0 / double(N_PHI_BINS);
    for (int i = 0; i < N_PHI_BINS; ++i) {
        v[i] = (i + 0.5) * step;
    }
    return v;
}

// Axis sets derived from xB/Q2/t ranges
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

// Per-bin maps in index space
struct PerBin {
    // key = (ix, iQ, it, ip)
    std::map<std::tuple<int,int,int,int>, double> lee_inb;
    std::map<std::tuple<int,int,int,int>, double> lee_out;
    std::map<std::tuple<int,int,int,int>, double> my_inb;
    std::map<std::tuple<int,int,int,int>, double> my_out;
};

static PerBin map_to_indices(const std::vector<BinRow>& rows,
                             const AxisSets& ax) {
    PerBin pb;

    for (const auto& r : rows) {
        const auto xb = std::make_pair(r.xBmin, r.xBmax);
        const int ix  = find_index(xb, ax.xB);
        if (ix < 0) continue;

        const auto& Q2s = ax.Q2_by_ix.at(ix);
        const auto& Ts  = ax.t_by_ix.at(ix);

        const int iQ = find_index({r.Q2min, r.Q2max}, Q2s);
        const int it = find_index({r.tmin,   r.tmax},  Ts);
        if (iQ < 0 || it < 0) continue;

        const int ip = phiBinFromDeg(r.phiavg);
        if (ip < 0) continue;

        auto key = std::make_tuple(ix, iQ, it, ip);

        pb.lee_inb[key] += r.lee_inb;
        pb.lee_out[key] += r.lee_out;
        pb.my_inb[key]  += r.my_inb;
        pb.my_out[key]  += r.my_out;
    }

    return pb;
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

static double compute_canvas_ymax(
    bool ratio_mode,
    const std::vector<std::pair<double,double>>& Q2s,
    const std::vector<std::pair<double,double>>& Ts,
    const std::function<bool(int,int,std::vector<double>&,std::vector<double>&)>& fillBoth
) {
    double ymax = 0.0;
    const int nrows = (int)Ts.size();
    const int ncols = (int)Q2s.size();

    for (int r = 0; r < nrows; ++r) {
        for (int ccol = 0; ccol < ncols; ++ccol) {
            std::vector<double> A(N_PHI_BINS, 0.0), B(N_PHI_BINS, 0.0);
            (void)fillBoth(ccol, r, A, B);

            if (!ratio_mode) {
                for (int ip = 0; ip < N_PHI_BINS; ++ip) {
                    ymax = std::max(ymax, A[ip]);
                    ymax = std::max(ymax, B[ip]);
                }
            } else {
                for (int ip = 0; ip < N_PHI_BINS; ++ip) {
                    const double NL = B[ip];
                    const double NH = A[ip];
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
                            const std::function<bool(int,int,std::vector<double>&,std::vector<double>&)>& fillBoth,
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

    const auto phiC = phiCentersDeg();
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

            std::vector<double> A(N_PHI_BINS, 0.0), B(N_PHI_BINS, 0.0);
            (void)fillBoth(ccol, r, A, B);

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

            if (draw_ratio_only) {
                std::vector<double> x, y, ey;
                x.reserve(N_PHI_BINS);
                y.reserve(N_PHI_BINS);
                ey.reserve(N_PHI_BINS);

                for (int ip = 0; ip < N_PHI_BINS; ++ip) {
                    const double NL = B[ip];
                    const double NH = A[ip];
                    if (NL <= 0.0) continue;
                    const double R = (NH <= 0.0) ? 0.0 : NH / NL;
                    double eR = 0.0;
                    if (NH > 0.0) {
                        eR = R * std::sqrt(1.0 / NH + 1.0 / NL);
                    }
                    x.push_back(phiC[ip]);
                    y.push_back(R);
                    ey.push_back(eR);
                }

                graph_pe1(x, y, ey, 20, black);

                TLine* one = new TLine(0.0, 1.0, 360.0, 1.0);
                one->SetLineStyle(2);
                one->SetLineWidth(2);
                one->SetLineColor(orange);
                one->Draw("SAME");
            } else {
                std::vector<double> eyA(N_PHI_BINS, 0.0), eyB(N_PHI_BINS, 0.0);
                for (int ip = 0; ip < N_PHI_BINS; ++ip) {
                    eyA[ip] = (A[ip] > 0.0) ? std::sqrt(A[ip]) : 0.0;
                    eyB[ip] = (B[ip] > 0.0) ? std::sqrt(B[ip]) : 0.0;
                }
                graph_pe1(phiC, A, eyA, 20, black);   // Hayward (pass-2)
                graph_pe1(phiC, B, eyB, 24, orange);  // Lee (pass-1)
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
    auto H = build_header_index(header);

    const std::vector<std::string> required = {
        "bin index", "valid bin",
        "xBmin", "xBmax",
        "Q2min", "Q2max",
        "tmin", "tmax",
        "phiavg"
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

        BinRow r;
        r.bin_index = ToInt(get_col_ref(cols, H, "bin index"));
        r.xBmin     = ToDouble(get_col_ref(cols, H, "xBmin"));
        r.xBmax     = ToDouble(get_col_ref(cols, H, "xBmax"));
        r.Q2min     = ToDouble(get_col_ref(cols, H, "Q2min"));
        r.Q2max     = ToDouble(get_col_ref(cols, H, "Q2max"));
        r.tmin      = ToDouble(get_col_ref(cols, H, "tmin"));
        r.tmax      = ToDouble(get_col_ref(cols, H, "tmax"));
        r.phiavg    = ToDouble(get_col_ref(cols, H, "phiavg"));

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

    // 3) Build axis sets and per-bin maps
    AxisSets ax = build_axes_from_rows(rows);
    PerBin pb   = map_to_indices(rows, ax);

    info("Axis xB bins: " + std::to_string(ax.xB.size()));

    auto make_fillBoth = [&](const std::map<std::tuple<int,int,int,int>, double>& ours,
                             const std::map<std::tuple<int,int,int,int>, double>& lee_side,
                             int ix) {
        return [&, ix](int iQcol, int irow,
                       std::vector<double>& A,
                       std::vector<double>& B)->bool {
            bool any = false;
            for (int ip = 0; ip < N_PHI_BINS; ++ip) {
                auto key = std::make_tuple(ix, iQcol, irow, ip);
                double a = 0.0;
                double b = 0.0;

                auto itA = ours.find(key);
                if (itA != ours.end()) a = itA->second;

                auto itB = lee_side.find(key);
                if (itB != lee_side.end()) b = itB->second;

                A[ip] = a;
                B[ip] = b;
                any  |= (a > 0.0 || b > 0.0);
            }
            return any;
        };
    };

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

            auto fillBoth_inb = make_fillBoth(pb.my_inb, pb.lee_inb, ix);

            const std::string f_counts =
                (fs::path(output_base_dir) / Form("raw_counts_fa18_inb_xB_%d.png", ix)).string();
            const std::string f_ratio  =
                (fs::path(output_base_dir) / Form("raw_ratio_fa18_inb_xB_%d.png",  ix)).string();

            draw_one_canvas(title_counts, Q2s, Ts, fillBoth_inb, f_counts, /*draw_ratio_only=*/false);
            draw_one_canvas(title_ratio,  Q2s, Ts, fillBoth_inb, f_ratio,   /*draw_ratio_only=*/true);

            info("Saved: " + f_counts);
            info("Saved: " + f_ratio);
        }

        // Fa18 Out
        {
            const std::string title_counts =
                Form("Raw yields: Fa18 Out   x_{B} #in [%.3g, %.3g]", xb_lo, xb_hi);
            const std::string title_ratio  =
                Form("Raw yields ratio (Hayward/Lee): Fa18 Out   x_{B} #in [%.3g, %.3g]", xb_lo, xb_hi);

            auto fillBoth_out = make_fillBoth(pb.my_out, pb.lee_out, ix);

            const std::string f_counts =
                (fs::path(output_base_dir) / Form("raw_counts_fa18_out_xB_%d.png", ix)).string();
            const std::string f_ratio  =
                (fs::path(output_base_dir) / Form("raw_ratio_fa18_out_xB_%d.png",  ix)).string();

            draw_one_canvas(title_counts, Q2s, Ts, fillBoth_out, f_counts, /*draw_ratio_only=*/false);
            draw_one_canvas(title_ratio,  Q2s, Ts, fillBoth_out, f_ratio,   /*draw_ratio_only=*/true);

            info("Saved: " + f_counts);
            info("Saved: " + f_ratio);
        }
    }
}