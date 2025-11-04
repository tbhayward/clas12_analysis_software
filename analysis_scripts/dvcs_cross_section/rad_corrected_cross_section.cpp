#include "rad_corrected_cross_section.h"

#include <TCanvas.h>
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
#include <vector>

#include <nlohmann/json.hpp>
using json = nlohmann::json;
namespace fs = std::filesystem;

namespace {
constexpr int    N_PHI_BINS = 12;
constexpr double TWO_PI     = 2.0 * M_PI;

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

static inline std::vector<double> phiCentersDeg() {
    std::vector<double> d(N_PHI_BINS);
    const double step = 360.0 / double(N_PHI_BINS);
    for (int i=0;i<N_PHI_BINS;++i) d[i] = (i+0.5)*step;
    return d;
}
static void drawDegreeTicks(double xmin, double ymin, double xmax, double labelSize){
    TGaxis* ax = new TGaxis(xmin, ymin, xmax, ymin, 0.0, 360.0, 4, "");
    ax->SetLabelFont(42);
    ax->SetLabelSize(labelSize);
    ax->SetLabelOffset(0.012);
    ax->SetTitle("");
    ax->SetTickSize(0.02);
    ax->Draw();
}

// bin helpers
static std::vector<std::pair<double,double>> uniqueRanges(const std::vector<Binning>& scheme, char which) {
    std::set<std::pair<double,double>> s;
    for (const auto& b : scheme) {
        if (which == 'x') s.emplace(b.xBmin, b.xBmax);
        else if (which == 'Q') s.emplace(b.Q2min, b.Q2max);
        else if (which == 't') s.emplace(b.tmin, b.tmax);
    }
    return std::vector<std::pair<double,double>>(s.begin(), s.end());
}
static inline int findIndex(const std::pair<double,double>& range,
                            const std::vector<std::pair<double,double>>& ranges) {
    for (int i=0;i<(int)ranges.size();++i) if (ranges[i]==range) return i;
    return -1;
}

// ---------- I/O helpers ----------
static json load_json(const std::string& path) {
    std::ifstream f(path);
    if (!f.is_open()) {
        std::cerr << "[json] Failed to open " << path << "\n";
        return json();
    }
    json j; f >> j; return j;
}

static inline bool has_bins12(const json& h){
    return h.contains("phi") && h.contains("xsec") && h.contains("xsec_err")
        && h["phi"].size()==N_PHI_BINS && h["xsec"].size()==N_PHI_BINS && h["xsec_err"].size()==N_PHI_BINS;
}
static inline bool has_rc12(const json& r){
    return r.contains("phi") && r.contains("rc")
        && r["phi"].size()==N_PHI_BINS && r["rc"].size()==N_PHI_BINS;
}

// safe read double from json index
static inline double jgetd(const json& a, size_t i, double def=0.0){
    try { return a[i].get<double>(); } catch(...) { return def; }
}

// Helper function to calculate log-scale y-axis range for a canvas
static std::pair<double, double> calculateYRangeForCanvas(const json& jout_bins, 
                                                         const json& j_unc_bins,
                                                         const std::vector<std::string>& bin_keys,
                                                         bool overlay) {
    double global_min = 1e10;
    double global_max = -1e10;
    
    for (const auto& bkey : bin_keys) {
        if (!jout_bins.contains(bkey)) continue;
        const json& jb_corr = jout_bins[bkey];
        
        // Check corrected data
        auto checkGraphRange = [&](const json& h) {
            if (!h.contains("xsec") || !h.contains("xsec_err")) return;
            const auto& yp = h["xsec"];
            const auto& ep = h["xsec_err"];
            for (int i = 0; i < N_PHI_BINS; ++i) {
                double y = jgetd(yp, i, 0.0);
                double e = jgetd(ep, i, 0.0);
                if (y > 0) {
                    global_min = std::min(global_min, std::max(1e-10, y - e));
                    global_max = std::max(global_max, y + e);
                }
            }
        };
        
        if (jb_corr.contains("helicity_plus")) 
            checkGraphRange(jb_corr["helicity_plus"]);
        if (jb_corr.contains("helicity_minus")) 
            checkGraphRange(jb_corr["helicity_minus"]);
        
        // Check uncorrected data if overlay
        if (overlay && j_unc_bins.contains(bkey)) {
            const json& jb_unc = j_unc_bins[bkey];
            if (jb_unc.contains("helicity_plus")) 
                checkGraphRange(jb_unc["helicity_plus"]);
            if (jb_unc.contains("helicity_minus")) 
                checkGraphRange(jb_unc["helicity_minus"]);
        }
    }
    
    // Set defaults if no valid data found
    if (global_max <= 0) {
        global_min = 1e-4;
        global_max = 1.0;
    } else {
        // Round down to nearest power of 10 below min, but set a floor of 1e-4
        double calculated_min = std::pow(10.0, std::floor(std::log10(global_min)));
        global_min = std::max(1e-4, calculated_min); // Use 1e-4 as minimum, or calculated if higher
        
        // Round up to nearest power of 10 above max
        global_max = std::pow(10.0, std::ceil(std::log10(global_max)));
        
        // Add some padding
        global_min *= 0.5;
        global_max *= 2.0;
        
        // Ensure we don't go below our global minimum after padding
        global_min = std::max(1e-4, global_min);
    }
    
    return {global_min, global_max};
}

} // anon

// ------------------------------------------------------------
// Main
// ------------------------------------------------------------
void compute_rad_corrected_cross_sections(
    const std::vector<Binning>& binning_scheme,
    const std::string& uncorrected_xsec_json_dir, // e.g. "output/jsons"
    const std::string& radcorr_json_dir,          // e.g. "output/jsons"
    const std::string& output_dir)                // e.g. "output/rad_corrected_cross_section"
{
    // dirs
    fs::create_directories(output_dir);
    fs::create_directories(fs::path(output_dir)/"plots");
    fs::create_directories(uncorrected_xsec_json_dir);
    fs::create_directories(radcorr_json_dir); // <<— JSONs will be saved here now

    // binning helpers for plotting layout
    const auto xB_bins = uniqueRanges(binning_scheme, 'x');
    const auto Q2_all  = uniqueRanges(binning_scheme, 'Q');
    const auto t_all   = uniqueRanges(binning_scheme, 't');

    // beam energies we know about; skip gracefully if files are missing
    const std::vector<std::string> energies = {"10.59","10.60","10.2"};
    const auto PHI_DEG = phiCentersDeg();

    for (const auto& E : energies) {
        const std::string un_path = (fs::path(uncorrected_xsec_json_dir) / ("uncorrected_xsec_" + E + ".json")).string();
        const std::string rc_path = (fs::path(radcorr_json_dir) / ("radiative_corrections_group_" + E + ".json")).string();

        json j_unc = load_json(un_path);
        json j_rc  = load_json(rc_path);

        if (j_unc.empty()) { std::cerr << "[radxsec] NOTE: missing uncorrected file for E="<<E<<" -> "<<un_path<<"\n"; continue; }
        if (j_rc.empty())  { std::cerr << "[radxsec] NOTE: missing RC file for E="<<E<<" -> "<<rc_path<<"\n"; continue; }

        if (!j_unc.contains("bins") || !j_rc.contains("bins")) {
            std::cerr << "[radxsec] Malformed JSON(s) for E="<<E<<" (no 'bins'). Skipping.\n";
            continue;
        }

        // output JSON container (corrected)
        json jout;
        jout["energy"] = E;
        jout["bins"]   = json::object();

        // Corrected = Uncorrected * (Born/Rad)
        for (auto it = j_unc["bins"].begin(); it != j_unc["bins"].end(); ++it) {
            const std::string bkey = it.key();
            const json& cell_unc = it.value();
            if (!j_rc["bins"].contains(bkey)) continue; // need matching RC
            const json& cell_rc = j_rc["bins"][bkey];

            if (!has_rc12(cell_rc)) continue;
            const auto& rc_phi = cell_rc["rc"];
            const bool has_rc_err = cell_rc.contains("rc_err") && cell_rc["rc_err"].size()==N_PHI_BINS;
            const auto& rc_err = has_rc_err ? cell_rc["rc_err"] : json::array();

            auto corrOneHel = [&](const char* node)->json{
                json out;
                if (!cell_unc.contains(node)) return out;
                const json& h = cell_unc[node];
                if (!has_bins12(h)) return out;

                const auto& ph = h["phi"];
                const auto& xs = h["xsec"];
                const auto& es = h["xsec_err"];

                std::vector<double> phi(N_PHI_BINS), y(N_PHI_BINS), e(N_PHI_BINS);
                for (int ip=0; ip<N_PHI_BINS; ++ip){
                    phi[ip] = jgetd(ph, ip, PHI_DEG[ip]);
                    const double x  = jgetd(xs, ip, 0.0);
                    const double ex = std::max(0.0, jgetd(es, ip, 0.0));
                    const double r  = std::max(0.0, jgetd(rc_phi, ip, 1.0));
                    const double er = has_rc_err ? std::max(0.0, jgetd(rc_err, ip, 0.0)) : 0.0;

                    // y' = r * x ;    sigma_y'^2 = (r*sigma_x)^2 + (x*sigma_r)^2
                    const double yc = r * x;
                    double ec = std::sqrt( (r*ex)*(r*ex) + (x*er)*(x*er) );
                    if (!std::isfinite(ec)) ec = 0.0;

                    y[ip] = yc; e[ip] = ec;
                }

                out["phi"]      = phi;
                out["xsec"]     = y;
                out["xsec_err"] = e;
                return out;
            };

            json hp = corrOneHel("helicity_plus");
            json hm = corrOneHel("helicity_minus");
            if (hp.is_null() && hm.is_null()) continue;

            json rc_used;
            rc_used["phi"]    = cell_rc["phi"];
            rc_used["rc"]     = rc_phi;
            if (has_rc_err) rc_used["rc_err"] = rc_err;

            jout["bins"][bkey] = {
                {"helicity_plus",  hp},
                {"helicity_minus", hm},
                {"rc_used",        rc_used}
            };
        }

        // Save corrected JSON — now in the universal output/jsons directory
        const std::string out_json = (fs::path(radcorr_json_dir)/("rad_corrected_xsec_"+E+".json")).string();
        {
            std::ofstream ofs(out_json);
            ofs << std::setw(2) << jout << "\n";
        }
        std::cout << "[radxsec] Wrote " << out_json << "\n";

        // ---------------- PLOTS (two variants) ----------------
        const auto Q2_all = uniqueRanges(binning_scheme, 'Q');
        const auto t_all  = uniqueRanges(binning_scheme, 't');

        auto drawHelFrom = [&](const json& jb, const char* node, int mstyle, int color){
            if (!jb.contains(node)) return (TGraphErrors*)nullptr;
            const auto& h = jb[node];
            if (!h.contains("phi") || !h.contains("xsec") || !h.contains("xsec_err")) return (TGraphErrors*)nullptr;
            const auto& xp = h["phi"];
            const auto& yp = h["xsec"];
            const auto& ep = h["xsec_err"];
            if (xp.size()!=N_PHI_BINS || yp.size()!=N_PHI_BINS || ep.size()!=N_PHI_BINS) return (TGraphErrors*)nullptr;
            std::vector<double> x(N_PHI_BINS), y(N_PHI_BINS), e(N_PHI_BINS);
            for (int i=0;i<N_PHI_BINS;++i){
                x[i]=xp[i].get<double>(); y[i]=yp[i].get<double>(); e[i]=std::max(1e-12, ep[i].get<double>());
            }
            auto* gr = new TGraphErrors(N_PHI_BINS, x.data(), y.data(), nullptr, e.data());
            gr->SetMarkerStyle(mstyle);
            gr->SetMarkerSize(1.0);
            gr->SetLineWidth(2);
            gr->SetLineColor(color);
            gr->SetMarkerColor(color);
            gr->Draw("P SAME");
            return gr;
        };

        auto makeCanvas = [&](bool overlay){
            // overlay=false: corrected-only
            // overlay=true : corrected + uncorrected together
            for (int ix = 0; ix < (int)xB_bins.size(); ++ix) {
                const auto xb = xB_bins[ix];

                // slice-specific lists
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
                cname << "c_radxsec_"<<E<<"_xB"<<ix<<(overlay ? "_overlay" : "_corr");
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
                head.SetTextSize(0.30);
                std::ostringstream tit;
                tit << (overlay ? "unc vs. rad corrected d#sigma/d#phi"
                                : "rad corrected d#sigma/d#phi");
                tit << Form("  %s GeV x_{B} #in [%.2g, %.2g]", E.c_str(), xb.first, xb.second);
                head.DrawLatex(0.5, 0.55, tit.str().c_str());

                // First pass: collect all bin keys for this xB slice to calculate y-range
                std::vector<std::string> canvas_bin_keys;
                for (int r=0; r<nrows; ++r) {
                    const int it_global = findIndex(t_slice[r], t_all);
                    if (it_global < 0) continue;
                    for (int cc=0; cc<ncols; ++cc) {
                        const int iQ_global = findIndex(Q2_slice[cc], Q2_all);
                        if (iQ_global < 0) continue;
                        const std::string bkey = "(" + std::to_string(ix) + "," + std::to_string(iQ_global) + "," + std::to_string(it_global) + ")";
                        canvas_bin_keys.push_back(bkey);
                    }
                }

                // Calculate y-range for this entire canvas
                auto [ymin_canvas, ymax_canvas] = calculateYRangeForCanvas(
                    jout["bins"], j_unc["bins"], canvas_bin_keys, overlay);

                // Second pass: draw each panel
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

                        // Use canvas-specific y-range
                        TH1* frame = gPad->DrawFrame(0.0, ymin_canvas, 360.0, ymax_canvas);
                        frame->GetXaxis()->SetLabelSize(0.0001);
                        frame->GetXaxis()->SetTitle("#phi (deg)");
                        frame->GetYaxis()->SetTitle(overlay ? "d#sigma/d#phi" : "d#sigma/d#phi (rad-corr.)");
                        frame->GetXaxis()->CenterTitle();
                        frame->GetYaxis()->CenterTitle();
                        frame->GetXaxis()->SetNdivisions(505);
                        frame->GetXaxis()->SetTitleSize(0.060);
                        frame->GetYaxis()->SetTitleSize(0.060);
                        frame->GetYaxis()->SetLabelSize(0.048);
                        frame->GetXaxis()->SetTitleOffset(1.25);
                        frame->GetYaxis()->SetTitleOffset(1.35);

                        drawDegreeTicks(gPad->GetUxmin(), gPad->GetUymin(), gPad->GetUxmax(), 0.048);

                        const std::string bkey = "(" + std::to_string(ix) + "," + std::to_string(iQ_global) + "," + std::to_string(it_global) + ")";
                        if (!jout["bins"].contains(bkey)) continue;
                        const json& jb_corr = jout["bins"][bkey];

                        // Optional uncorrected overlay
                        const bool have_unc = overlay && j_unc["bins"].contains(bkey);
                        const json& jb_unc = have_unc ? j_unc["bins"][bkey] : json::object();

                        // draw order: uncorrected first (grey), then corrected (color)
                        TGraphErrors* gup = nullptr;
                        TGraphErrors* gum = nullptr;
                        if (have_unc) {
                            gup = drawHelFrom(jb_unc,  "helicity_plus",  24 /*open circle*/, kGray+2);
                            gum = drawHelFrom(jb_unc,  "helicity_minus", 26 /*open triangle*/, kGray+2);
                        }
                        TGraphErrors* gcp = drawHelFrom(jb_corr, "helicity_plus",  20 /*filled circle*/, kBlue+1);
                        TGraphErrors* gcm = drawHelFrom(jb_corr, "helicity_minus", 25 /*filled square*/, kRed+1);

                        TLatex lab;
                        lab.SetNDC(); lab.SetTextSize(0.040); lab.SetTextAlign(11);
                        lab.SetTextFont(48);
                        lab.DrawLatex(0.15, 0.94,
                            Form("Q^{2} #in [%.2g, %.2g], -t #in [%.2g, %.2g]",
                                 Q2_slice[cc].first, Q2_slice[cc].second,
                                 t_slice[r].first,   t_slice[r].second));

                        if (gcp || gcm || gup || gum){
                            double x1=0.44, y1=0.68, x2=0.90, y2=0.93;
                            TLegend* leg = new TLegend(x1,y1,x2,y2);
                            leg->SetBorderSize(1);
                            leg->SetLineColor(kBlack);
                            leg->SetFillColor(kWhite);
                            leg->SetFillStyle(1001);
                            leg->SetTextFont(42);
                            leg->SetTextSize(0.042);
                            if (gcp) leg->AddEntry(gcp, "corrected  + helicity", "lep");
                            if (gcm) leg->AddEntry(gcm, "corrected  - helicity", "lep");
                            if (gup) leg->AddEntry(gup, "uncorrected + helicity", "lep");
                            if (gum) leg->AddEntry(gum, "uncorrected - helicity", "lep");
                            leg->Draw();
                        }
                    }
                }

                const std::string outP =
                    (fs::path(output_dir)/"plots"/(std::string("rad_corrected_xsec_")+E+"_xB_"+std::to_string(ix)+(overlay ? "_overlay.png" : ".png"))).string();
                c->SaveAs(outP.c_str());
                delete c;
            }
        };

        // Save both variants
        makeCanvas(false); // corrected-only
        makeCanvas(true);  // overlay
    }

    std::cout << "[radxsec] Finished radiatively corrected cross-section generation.\n";
}