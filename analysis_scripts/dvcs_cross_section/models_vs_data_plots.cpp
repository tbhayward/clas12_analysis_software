#include "models_vs_data_plots.h"

#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TPad.h>
#include <TH1.h>
#include <TAxis.h>
#include <TGaxis.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

// OpenMP only for model sampling; ROOT drawing remains single-threaded.
#ifdef _OPENMP
#include <omp.h>
#endif

#include <nlohmann/json.hpp>
using json = nlohmann::json;
namespace fs = std::filesystem;

namespace {
constexpr int    N_PHI_BINS = 12;

// ------------ style ------------
struct StyleInit {
    StyleInit() {
        gStyle->SetOptTitle(0);
        gStyle->SetOptStat(0);
        gStyle->SetLineWidth(2);
        gStyle->SetFrameLineWidth(2);
        gStyle->SetPadTickX(1);
        gStyle->SetPadTickY(1);
        const int rf = 42;
        gStyle->SetTitleFont(rf, "XYZ");
        gStyle->SetLabelFont(rf, "XYZ");
        gStyle->SetTextFont(rf);
    }
} _style_guard;

// Helpers from your plotting conventions
static inline std::vector<double> phiCentersDeg() {
    std::vector<double> d(N_PHI_BINS);
    const double step = 360.0 / double(N_PHI_BINS);
    for (int i=0;i<N_PHI_BINS;++i) d[i] = (i+0.5)*step;
    return d;
} //enddef

static void drawDegreeTicks(double xmin, double ymin, double xmax, double labelSize) {
    TGaxis* ax = new TGaxis(xmin, ymin, xmax, ymin, 0.0, 360.0, 4, "");
    ax->SetLabelFont(42);
    ax->SetLabelSize(labelSize);
    ax->SetLabelOffset(0.012);
    ax->SetTitle("");
    ax->SetTickSize(0.02);
    ax->Draw();
} //enddef

// bin helpers (duplicate of your other modules to stay standalone)
static std::vector<std::pair<double,double>> uniqueRanges(const std::vector<Binning>& scheme, char which) {
    std::set<std::pair<double,double>> s;
    for (const auto& b : scheme) {
        if      (which == 'x') s.emplace(b.xBmin, b.xBmax);
        else if (which == 'Q') s.emplace(b.Q2min, b.Q2max);
        else if (which == 't') s.emplace(b.tmin, b.tmax);
    } //endfor
    return std::vector<std::pair<double,double>>(s.begin(), s.end());
} //enddef

static inline int findIndex(const std::pair<double,double>& range,
                            const std::vector<std::pair<double,double>>& ranges) {
    for (int i=0;i<(int)ranges.size();++i) if (ranges[i]==range) return i; //endfor
    return -1;
} //enddef

// JSON helpers
static json load_json(const std::string& path) {
    std::ifstream f(path);
    if (!f.is_open()) {
        std::cerr << "[models] Failed to open " << path << "\n";
        return json();
    } //endif
    json j; f >> j; return j;
} //enddef

static inline bool has_bins12(const json& h){
    return h.contains("phi") && h.contains("xsec") && h.contains("xsec_err")
        && h["phi"].size()==N_PHI_BINS && h["xsec"].size()==N_PHI_BINS && h["xsec_err"].size()==N_PHI_BINS;
} //enddef

// safe read double from json index
static inline double jgetd(const json& a, size_t i, double def=0.0){
    try { return a[i].get<double>(); } catch(...) { return def; }
} //enddef

static inline double midpoint(double a, double b){ return 0.5*(a+b); } //enddef
static inline bool positive(double v){ return std::isfinite(v) && v>0.0; } //enddef

// Draw the data points (bin-centered cross sections)
static TGraphErrors* draw_data_hel(const json& h, int mstyle, int color, TH1* frame, double& ymax_accum) {
    if (!has_bins12(h)) return nullptr; //endif
    const auto& xp = h["phi"];
    const auto& yp = h["xsec"];
    const auto& ep = h["xsec_err"];
    std::vector<double> x(N_PHI_BINS), y(N_PHI_BINS), e(N_PHI_BINS);
    for (int i=0;i<N_PHI_BINS;++i){
        x[i]=xp[i].get<double>();
        y[i]=yp[i].get<double>();
        e[i]=std::max(1e-12, ep[i].get<double>());
        ymax_accum = std::max(ymax_accum, y[i]+e[i]);
    } //endfor
    auto* gr = new TGraphErrors(N_PHI_BINS, x.data(), y.data(), nullptr, e.data());
    gr->SetMarkerStyle(mstyle);
    gr->SetMarkerSize(1.0);
    gr->SetLineWidth(2);
    gr->SetLineColor(color);
    gr->SetMarkerColor(color);
    gr->Draw("P SAME");
    return gr;
} //enddef

// Sample a model curve densely in phi (0..360) at bin centers (xB,Q2,t)
struct ModelCurves {
    std::vector<double> phi_deg;  // 0..360
    std::vector<double> y_plus;   // + helicity
    std::vector<double> y_minus;  // - helicity
    std::vector<double> y_bh;     // BH (helicity-independent)
};

static ModelCurves sample_models_dense(
    int phi_dense,
    double Ebeam,
    double xBmin, double xBmax,
    double Q2min, double Q2max,
    double tmin_pos, double tmax_pos,
    const ModelPaths& paths,
    bool vgg_globalfit
) {
    ModelCurves mc;
    if (phi_dense < 2) phi_dense = 2;
    mc.phi_deg.resize(phi_dense);
    mc.y_plus.resize(phi_dense);
    mc.y_minus.resize(phi_dense);
    mc.y_bh.resize(phi_dense);

    const double xB_c   = midpoint(xBmin, xBmax);
    const double Q2_c   = midpoint(Q2min, Q2max);
    const double tpos_c = midpoint(tmin_pos, tmax_pos);

    // Hard cap: use at most 5 threads here, no matter the machine.
    #ifdef _OPENMP
    const int user_cap = 5;
    int sys_cap = omp_get_max_threads();
    int nthreads = std::min(user_cap, std::max(1, sys_cap));
    omp_set_num_threads(nthreads);
    #endif

    // Sample 0..360 inclusive (wrap last point to 360 exactly)
    #pragma omp parallel for if(phi_dense > 72) schedule(static)
    for (int i=0; i<phi_dense; ++i) {
        double ph = (phi_dense == 1) ? 0.0 : (360.0 * double(i) / double(phi_dense-1));
        mc.phi_deg[i] = ph;

        // KM15, VGG, BH are combined later by line style; here we store final choices:
        // Convention: we will compute VGG and KM15 in the draw step; to get a y-range that
        // accounts for BH too, compute BH here as well.
        // If your BH function name differs, adjust here.
        double ybh = 0.0;
        try {
            ybh = bh_xs(xB_c, Q2_c, tpos_c, ph, Ebeam, Helicity::Plus, paths); // helicity-independent
        } catch (...) {
            ybh = 0.0;
        }
        mc.y_bh[i] = positive(ybh) ? ybh : 0.0;

        // Placeholders; will be overwritten in draw with actual VGG/KM15 samples if desired.
        mc.y_plus[i]  = 0.0;
        mc.y_minus[i] = 0.0;
    } //endfor

    return mc;
} //enddef

// Compute a specific model at (phi, +) and (phi, -) for a given label
enum class ModelTag { VGG, KM15 };

static inline double eval_model(ModelTag tag, double xB, double Q2, double tpos, double ph_deg,
                                double Ebeam, Helicity hel, const ModelPaths& paths, bool vgg_globalfit) {
    if (tag == ModelTag::KM15) {
        return km15_xs(xB, Q2, tpos, ph_deg, Ebeam, hel, paths);
    } else {
        return vgg_xs(xB, Q2, tpos, ph_deg, Ebeam, hel, paths, vgg_globalfit);
    }
} //enddef

// Compute canvas-wide y-range including data points and model envelopes
static std::pair<double,double> computeYRangeCanvas(
    const json& jout_bins,                                       // data (bin-centered)
    const std::vector<std::string>& bin_keys,
    const std::map<std::string, std::pair<double,double>>& model_envelopes  // per-panel (min,max) from model sampling
) {
    double global_min = 1e10;
    double global_max = -1e10;

    auto collectData = [&](const json& h) {
        if (!h.contains("xsec") || !h.contains("xsec_err")) return;
        const auto& yp = h["xsec"];
        const auto& ep = h["xsec_err"];
        for (int i=0;i<N_PHI_BINS;++i) {
            double y = jgetd(yp, i, 0.0);
            double e = jgetd(ep, i, 0.0);
            if (y > 0.0) {
                global_min = std::min(global_min, std::max(1e-12, y - e));
                global_max = std::max(global_max, y + e);
            }
        }
    };

    for (const auto& bkey : bin_keys) {
        if (!jout_bins.contains(bkey)) continue;
        const json& jb = jout_bins[bkey];
        if (jb.contains("helicity_plus"))  collectData(jb["helicity_plus"]);
        if (jb.contains("helicity_minus")) collectData(jb["helicity_minus"]);

        auto it = model_envelopes.find(bkey);
        if (it != model_envelopes.end()) {
            global_min = std::min(global_min, it->second.first);
            global_max = std::max(global_max, it->second.second);
        }
    }

    // Defaults if nothing valid
    if (!(global_max > 0.0)) return {1e-4, 1.0};

    // Pad and round
    double ymin = std::pow(10.0, std::floor(std::log10(std::max(1e-12, global_min))));
    ymin = std::max(1e-4, ymin*0.5);
    double ymax = std::pow(10.0, std::ceil(std::log10(global_max)));
    ymax *= 2.0;

    return {ymin, ymax};
} //enddef

} // anon

void plot_models_vs_bincentered(
    const std::vector<Binning>& binning_scheme,
    const std::string& bincenter_json_dir,
    const std::string& output_dir,
    const ModelPaths& paths,
    bool vgg_globalfit,
    int phi_dense
) {
    fs::create_directories(output_dir);
    fs::create_directories(fs::path(output_dir)/"plots");

    const auto xB_bins = uniqueRanges(binning_scheme, 'x');
    const auto Q2_all  = uniqueRanges(binning_scheme, 'Q');
    const auto t_all   = uniqueRanges(binning_scheme, 't');

    const std::vector<std::string> energies = {"10.59","10.60","10.2"};
    const auto PHI_DEG = phiCentersDeg();

    for (const auto& E : energies) {
        const std::string bc_path = (fs::path(bincenter_json_dir) / ("bin_centered_xsec_" + E + ".json")).string();
        json j_bc = load_json(bc_path);
        if (j_bc.empty() || !j_bc.contains("bins")) {
            std::cerr << "[models] NOTE: missing or malformed bin-centered file for E=" << E
                      << " -> " << bc_path << " (skipping)\n";
            continue;
        }

        // For each xB slice, build a canvas of Q2 x t panels
        for (int ix = 0; ix < (int)xB_bins.size(); ++ix) {
            const auto xb = xB_bins[ix];

            // Discover which Q2,t actually appear in this xB slice from the scheme
            std::set<std::pair<double,double>> qs, ts;
            for (const auto& b : binning_scheme) {
                if (std::make_pair(b.xBmin,b.xBmax) == xb) {
                    qs.emplace(b.Q2min,b.Q2max);
                    ts.emplace(b.tmin,b.tmax);
                }
            }
            std::vector<std::pair<double,double>> Q2_slice(qs.begin(), qs.end());
            std::vector<std::pair<double,double>> t_slice(ts.begin(), ts.end());
            if (Q2_slice.empty() || t_slice.empty()) continue;

            const int nrows = (int)t_slice.size();
            const int ncols = (int)Q2_slice.size();
            const int W = 280*ncols + 160;
            const int H = 240*nrows + 170;

            std::ostringstream cname;
            cname << "c_models_"<<E<<"_xB"<<ix;
            TCanvas* c = new TCanvas(cname.str().c_str(), cname.str().c_str(), W, H);

            TPad* pTop  = new TPad("pTop","pTop", 0.0, 0.915, 1.0, 1.0);
            pTop->SetFillStyle(0); pTop->SetBorderSize(0); pTop->Draw();

            TPad* pGrid = new TPad("pGrid","pGrid", 0.0, 0.00, 1.0, 0.915);
            pGrid->SetFillStyle(0); pGrid->SetBorderSize(0); pGrid->Draw();
            pGrid->cd();
            pGrid->Divide(ncols, nrows, 0.0001, 0.0001);

            // Title
            pTop->cd();
            TLatex head;
            head.SetNDC(); head.SetTextAlign(22);
            head.SetTextFont(42);
            head.SetTextSize(0.22);
            std::ostringstream tit;
            tit << "models vs. bin-centered d#sigma/d#phi"
                << Form("  %s GeV x_{B} #in [%.2g, %.2g]", E.c_str(), xb.first, xb.second);
            head.DrawLatex(0.5, 0.55, tit.str().c_str());

            // First pass over panels: pre-sample BH and collect model envelopes for y-range
            std::map<std::string, std::pair<double,double>> model_envelopes; // bkey -> (min,max)
            for (int r=0; r<nrows; ++r) {
                const int it_global = findIndex(t_slice[r], t_all);
                if (it_global < 0) continue;
                for (int cc=0; cc<ncols; ++cc) {
                    const int iQ_global = findIndex(Q2_slice[cc], Q2_all);
                    if (iQ_global < 0) continue;

                    const std::string bkey =
                        "(" + std::to_string(ix) + "," + std::to_string(iQ_global) + "," + std::to_string(it_global) + ")";

                    // Pre-sample BH for envelope (KM15/VGG envelopes are added during draw)
                    ModelCurves mc_bh = sample_models_dense(
                        phi_dense, std::stod(E),
                        xb.first, xb.second,
                        Q2_slice[cc].first, Q2_slice[cc].second,
                        t_slice[r].first,   t_slice[r].second,
                        paths, vgg_globalfit
                    );
                    double mmin = 1e10, mmax = -1e10;
                    for (size_t i=0;i<mc_bh.y_bh.size();++i) {
                        double v = mc_bh.y_bh[i];
                        if (positive(v)) {
                            mmin = std::min(mmin, v);
                            mmax = std::max(mmax, v);
                        }
                    }
                    if (mmax <= 0.0) { mmin = 1e10; mmax = -1e10; } // if BH had no data
                    model_envelopes[bkey] = {mmin, mmax};
                }
            }

            // Second pass: compute canvas y-range from data and BH envelopes
            std::vector<std::string> canvas_bin_keys;
            for (int r=0; r<nrows; ++r) {
                const int it_global = findIndex(t_slice[r], t_all);
                if (it_global < 0) continue;
                for (int cc=0; cc<ncols; ++cc) {
                    const int iQ_global = findIndex(Q2_slice[cc], Q2_all);
                    if (iQ_global < 0) continue;
                    canvas_bin_keys.emplace_back(
                        "(" + std::to_string(ix) + "," + std::to_string(iQ_global) + "," + std::to_string(it_global) + ")"
                    );
                }
            }
            auto [ymin_canvas, ymax_canvas] = computeYRangeCanvas(j_bc["bins"], canvas_bin_keys, model_envelopes);

            // Third pass: draw each panel
            for (int r=0; r<nrows; ++r) {
                const int it_global = findIndex(t_slice[r], t_all);
                if (it_global < 0) continue;
                for (int cc=0; cc<ncols; ++cc) {
                    const int iQ_global = findIndex(Q2_slice[cc], Q2_all);
                    if (iQ_global < 0) continue;

                    pGrid->cd(r*ncols + cc + 1);
                    gPad->SetGrid(1,1);
                    gPad->SetTopMargin(0.08);
                    gPad->SetBottomMargin(0.18);
                    gPad->SetLeftMargin(0.15);
                    gPad->SetRightMargin(0.10);
                    gPad->SetLogy();

                    TH1* frame = gPad->DrawFrame(0.0, ymin_canvas, 360.0, ymax_canvas);
                    frame->GetXaxis()->SetLabelSize(0.0001);
                    frame->GetXaxis()->SetTitle("#phi (deg)");
                    frame->GetYaxis()->SetTitle("d#sigma/d#phi (bin-centered)");
                    frame->GetXaxis()->CenterTitle();
                    frame->GetYaxis()->CenterTitle();
                    frame->GetXaxis()->SetNdivisions(505);
                    frame->GetXaxis()->SetTitleSize(0.060);
                    frame->GetYaxis()->SetTitleSize(0.060);
                    frame->GetYaxis()->SetLabelSize(0.048);
                    frame->GetXaxis()->SetTitleOffset(1.25);
                    frame->GetYaxis()->SetTitleOffset(1.35);

                    drawDegreeTicks(gPad->GetUxmin(), gPad->GetUymin(), gPad->GetUxmax(), 0.048);

                    const std::string bkey =
                        "(" + std::to_string(ix) + "," + std::to_string(iQ_global) + "," + std::to_string(it_global) + ")";
                    if (!j_bc["bins"].contains(bkey)) continue;
                    const json& jb = j_bc["bins"][bkey];

                    // 1) Data points (bin-centered)
                    double ymax_seen = 0.0;
                    TGraphErrors* gdp = nullptr;
                    TGraphErrors* gdm = nullptr;
                    if (jb.contains("helicity_plus"))
                        gdp = draw_data_hel(jb["helicity_plus"],  20 /*filled circle*/, kBlue+1, frame, ymax_seen);
                    if (jb.contains("helicity_minus"))
                        gdm = draw_data_hel(jb["helicity_minus"], 25 /*filled square*/, kRed+1,  frame, ymax_seen);

                    // 2) Model curves (KM15 dashed, VGG dotted), and BH dashed black.
                    //    Sample densely now (at bin centers).
                    const double xB_c   = midpoint(xb.first, xb.second);
                    const double Q2_c   = midpoint(Q2_slice[cc].first, Q2_slice[cc].second);
                    const double tpos_c = midpoint(t_slice[r].first,   t_slice[r].second);

                    // Dense phi vector
                    std::vector<double> ph(phi_dense);
                    std::vector<double> y_km_plus(phi_dense),  y_km_minus(phi_dense);
                    std::vector<double> y_vg_plus(phi_dense),  y_vg_minus(phi_dense);
                    std::vector<double> y_bh(phi_dense);

                    for (int i=0; i<phi_dense; ++i) {
                        ph[i] = (phi_dense == 1) ? 0.0 : (360.0 * double(i) / double(phi_dense-1));
                    }

                    // Compute (no drawing) — allow limited parallel sampling
                    #ifdef _OPENMP
                    omp_set_num_threads(std::min(5, std::max(1, omp_get_max_threads())));
                    #endif
                    #pragma omp parallel for if(phi_dense > 72) schedule(static)
                    for (int i=0; i<phi_dense; ++i) {
                        double p = ph[i];
                        // KM15
                        double km_p = km15_xs(xB_c, Q2_c, tpos_c, p, std::stod(E), Helicity::Plus,  paths);
                        double km_m = km15_xs(xB_c, Q2_c, tpos_c, p, std::stod(E), Helicity::Minus, paths);
                        y_km_plus[i]  = positive(km_p) ? km_p : 0.0;
                        y_km_minus[i] = positive(km_m) ? km_m : 0.0;
                        // VGG
                        double vg_p = vgg_xs(xB_c, Q2_c, tpos_c, p, std::stod(E), Helicity::Plus,  paths, vgg_globalfit);
                        double vg_m = vgg_xs(xB_c, Q2_c, tpos_c, p, std::stod(E), Helicity::Minus, paths, vgg_globalfit);
                        y_vg_plus[i]  = positive(vg_p) ? vg_p : 0.0;
                        y_vg_minus[i] = positive(vg_m) ? vg_m : 0.0;
                        // BH (helicity independent)
                        double yb = 0.0;
                        try { yb = bh_xs(xB_c, Q2_c, tpos_c, p, std::stod(E), Helicity::Plus, paths); }
                        catch (...) { yb = 0.0; }
                        y_bh[i] = positive(yb) ? yb : 0.0;
                    }

                    // Build TGraphs and draw (lines only)
                    auto makeLine = [&](const std::vector<double>& xv, const std::vector<double>& yv,
                                        int color, int lstyle)->TGraph* {
                        auto* gr = new TGraph((int)xv.size(), const_cast<double*>(xv.data()),
                                              const_cast<double*>(yv.data()));
                        gr->SetLineColor(color);
                        gr->SetLineWidth(3);
                        gr->SetLineStyle(lstyle);
                        gr->Draw("L SAME");
                        return gr;
                    };

                    // KM15 dashed: + blue, − red
                    TGraph* g_km_p = makeLine(ph, y_km_plus,  kBlue+1, 2);
                    TGraph* g_km_m = makeLine(ph, y_km_minus, kRed+1,  2);

                    // VGG dotted: + blue, − red
                    TGraph* g_vg_p = makeLine(ph, y_vg_plus,  kBlue+1, 3);
                    TGraph* g_vg_m = makeLine(ph, y_vg_minus, kRed+1,  3);

                    // BH dashed black
                    TGraph* g_bh   = makeLine(ph, y_bh, kBlack, 2);

                    // Panel label
                    TLatex lab;
                    lab.SetNDC(); lab.SetTextSize(0.07); lab.SetTextAlign(11);
                    lab.SetTextFont(42);
                    lab.DrawLatex(0.15, 0.94,
                        Form("Q^{2} #in [%.2g, %.2g], -t #in [%.2g, %.2g]",
                             Q2_slice[cc].first, Q2_slice[cc].second,
                             t_slice[r].first,   t_slice[r].second));

                    // Legend (bottom-left) — width per your tweak
                    if (gdp || gdm || g_km_p || g_km_m || g_vg_p || g_vg_m || g_bh) {
                        double x1=0.16, y1=0.18, x2=0.68, y2=0.42; // your updated width
                        TLegend* leg = new TLegend(x1,y1,x2,y2);
                        leg->SetBorderSize(1);
                        leg->SetLineColor(kBlack);
                        leg->SetFillColor(kWhite);
                        leg->SetFillStyle(1001);
                        leg->SetTextFont(42);
                        leg->SetTextSize(0.048);

                        if (gdp)    leg->AddEntry(gdp,   "+ helicity (data)", "lep");
                        if (gdm)    leg->AddEntry(gdm,   "- helicity (data)", "lep");
                        if (g_km_p) leg->AddEntry(g_km_p,"KM15 + (dashed)",   "l");
                        if (g_km_m) leg->AddEntry(g_km_m,"KM15 - (dashed)",   "l");
                        if (g_vg_p) leg->AddEntry(g_vg_p,"VGG + (dotted)",    "l");
                        if (g_vg_m) leg->AddEntry(g_vg_m,"VGG - (dotted)",    "l");
                        if (g_bh)   leg->AddEntry(g_bh,  "BH exact (dashed)", "l");
                        leg->Draw();
                    }
                }
            }

            const std::string outP =
                (fs::path(output_dir)/"plots"/(std::string("models_vs_bincentered_")+E+"_xB_"+std::to_string(ix)+".png")).string();
            c->SaveAs(outP.c_str());
            delete c;
        }
    }

    std::cout << "[models] Finished model vs. bin-centered plotting.\n";
}