#include "bin_centering_corrections.h"

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
} //enddef

static inline std::pair<double,double> phiEdgesDegForBin(int ip) {
    const double step = 360.0 / double(N_PHI_BINS);
    const double lo = ip * step;
    const double hi = (ip+1) * step;
    return {lo, hi};
} //enddef

static void drawDegreeTicks(double xmin, double ymin, double xmax, double labelSize){
    TGaxis* ax = new TGaxis(xmin, ymin, xmax, ymin, 0.0, 360.0, 4, "");
    ax->SetLabelFont(42);
    ax->SetLabelSize(labelSize);
    ax->SetLabelOffset(0.012);
    ax->SetTitle("");
    ax->SetTickSize(0.02);
    ax->Draw();
} //enddef

// bin helpers
static std::vector<std::pair<double,double>> uniqueRanges(const std::vector<Binning>& scheme, char which) {
    std::set<std::pair<double,double>> s;
    for (const auto& b : scheme) {
        if (which == 'x') s.emplace(b.xBmin, b.xBmax);
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

// ---------- I/O helpers ----------
static json load_json(const std::string& path) {
    std::ifstream f(path);
    if (!f.is_open()) {
        std::cerr << "[json] Failed to open " << path << "\n";
        return json();
    } //endif
    json j; f >> j; return j;
} //enddef

static void save_json(const json& j, const std::string& path) {
    std::ofstream ofs(path);
    if (!ofs.is_open()) {
        std::cerr << "[json] Failed to write " << path << "\n";
        return;
    } //endif
    ofs << std::setw(2) << j << "\n";
} //enddef

static inline bool has_bins12(const json& h){
    return h.contains("phi") && h.contains("xsec") && h.contains("xsec_err")
        && h["phi"].size()==N_PHI_BINS && h["xsec"].size()==N_PHI_BINS && h["xsec_err"].size()==N_PHI_BINS;
} //enddef

static inline bool has_bc_cell(const json& cell){
    return cell.contains("helicity_plus") || cell.contains("helicity_minus");
} //enddef

// safe read double from json index
static inline double jgetd(const json& a, size_t i, double def=0.0){
    try { return a[i].get<double>(); } catch(...) { return def; }
} //enddef

static inline double midpoint(double a, double b){ return 0.5*(a+b); } //enddef

static inline bool finite_pos(double v){ return std::isfinite(v) && v>0.0; } //enddef

// ------------------------------------------------------------
// Optimized helper: compute fbin for one (xB,Q2,t,phi-bin,helicity)
// ------------------------------------------------------------
static void compute_fbin_for_point(
    // bin ranges and centers
    double xBmin, double xBmax,
    double Q2min, double Q2max,
    double tmin_pos, double tmax_pos,   // positive convention for |t|
    double phi_center_deg,
    std::pair<double,double> phi_edges_deg,
    // model config
    int n_steps,
    double Ebeam,
    Helicity hel,
    const ModelPaths& paths,
    bool vgg_globalfit,
    ModelChoice model_choice,
    // outputs
    double& fbin_km15, double& fbin_vgg
) {
    fbin_km15 = 0.0;
    fbin_vgg  = 0.0;

    // centers
    const double xB_c   = midpoint(xBmin, xBmax);
    const double Q2_c   = midpoint(Q2min, Q2max);
    const double tpos_c = midpoint(tmin_pos, tmax_pos);

    // Only compute center values for models we're actually using
    double km15_center = 0.0;
    double vgg_center = 0.0;
    
    if (model_choice == ModelChoice::Both || model_choice == ModelChoice::KM15Only) {
        km15_center = km15_xs(xB_c, Q2_c, tpos_c, phi_center_deg, Ebeam, hel, paths);
    }
    if (model_choice == ModelChoice::Both || model_choice == ModelChoice::VGGOnly) {
        vgg_center  = vgg_xs(xB_c, Q2_c, tpos_c, phi_center_deg, Ebeam, hel, paths, vgg_globalfit);
    }

    // Early exit if we don't need to compute averages (no valid center values)
    bool need_km15 = (model_choice == ModelChoice::Both || model_choice == ModelChoice::KM15Only) && finite_pos(km15_center);
    bool need_vgg = (model_choice == ModelChoice::Both || model_choice == ModelChoice::VGGOnly) && finite_pos(vgg_center);
    
    if (!need_km15 && !need_vgg) {
        return;
    }

    // Use adaptive sub-binning: fewer steps for smaller bins or when n_steps is large
    const double xB_width = xBmax - xBmin;
    const double Q2_width = Q2max - Q2min;
    const double t_width = tmax_pos - tmin_pos;
    const double phi_width = phi_edges_deg.second - phi_edges_deg.first;
    
    // Adaptive step sizing: use fewer steps when the bin is small relative to typical ranges
    const int nx = (xB_width < 0.1) ? std::max(2, n_steps/2) : std::max(2, n_steps);
    const int nQ = (Q2_width < 1.0) ? std::max(2, n_steps/2) : std::max(2, n_steps);
    const int nt = (t_width < 0.1) ? std::max(2, n_steps/2) : std::max(2, n_steps);
    const int np = std::max(2, n_steps); // phi usually needs more points due to oscillations

    auto linspace = [](double lo, double hi, int n){
        std::vector<double> v(n);
        if (n==1) { v[0]=0.5*(lo+hi); return v; } //endif
        const double step = (hi - lo) / double(n-1);
        for (int i=0;i<n;++i) v[i] = lo + i*step; //endfor
        return v;
    }; //enddef

    const auto xs = linspace(xBmin, xBmax, nx);
    const auto Qs = linspace(Q2min, Q2max, nQ);
    const auto ts = linspace(tmin_pos, tmax_pos, nt);
    const auto ps = linspace(phi_edges_deg.first, phi_edges_deg.second, np);

    // accumulate over sub-bins
    double sum_km15 = 0.0; int cnt_km15 = 0;
    double sum_vgg  = 0.0; int cnt_vgg  = 0;

    int total_sub_bins = nx * nQ * nt * np;
    
    // Only show progress for very large computations
    if (total_sub_bins > 1000) {
        std::cout << "[bincenter]       Computing " << total_sub_bins << " sub-bins" << std::endl;
    }

    // Precompute values to avoid repeated function calls for the same parameters
    for (double xb : xs){
        for (double q2 : Qs){
            for (double tp : ts){
                // Compute both models at once for the same (xb, q2, tp) point across all phi values
                for (double ph : ps){
                    if (need_km15) {
                        const double vk = km15_xs(xb, q2, tp, ph, Ebeam, hel, paths);
                        if (finite_pos(vk)) { sum_km15 += vk; ++cnt_km15; } //endif
                    }
                    if (need_vgg) {
                        const double vv = vgg_xs(xb, q2, tp, ph, Ebeam, hel, paths, vgg_globalfit);
                        if (finite_pos(vv)) { sum_vgg  += vv; ++cnt_vgg; } //endif
                    }
                } //endfor
            } //endfor
        } //endfor
    } //endfor

    // averages
    if (need_km15 && cnt_km15 > 0) {
        const double avg_km15 = sum_km15 / double(cnt_km15);
        if (finite_pos(avg_km15)) fbin_km15 = km15_center / avg_km15;
    }
    
    if (need_vgg && cnt_vgg > 0) {
        const double avg_vgg = sum_vgg / double(cnt_vgg);
        if (finite_pos(avg_vgg)) fbin_vgg = vgg_center / avg_vgg;
    }
} //enddef

// ------------------------------------------------------------
// Plot helpers
// ------------------------------------------------------------
static TGraphErrors* draw_hel_curve(const json& h, int mstyle, int color, TH1* frame){
    if (!has_bins12(h)) return nullptr; //endif
    const auto& xp = h["phi"];
    const auto& yp = h["xsec"];
    const auto& ep = h["xsec_err"];
    std::vector<double> x(N_PHI_BINS), y(N_PHI_BINS), e(N_PHI_BINS);
    double ymax=0.0;
    for (int i=0;i<N_PHI_BINS;++i){
        x[i]=xp[i].get<double>();
        y[i]=yp[i].get<double>();
        e[i]=std::max(1e-12, ep[i].get<double>());
        ymax = std::max(ymax, y[i]+e[i]);
    } //endfor
    if (ymax > 0.0) frame->GetYaxis()->SetRangeUser(1e-4, std::max(1.0, ymax*1.5)); //endif
    auto* gr = new TGraphErrors(N_PHI_BINS, x.data(), y.data(), nullptr, e.data());
    gr->SetMarkerStyle(mstyle);
    gr->SetMarkerSize(1.0);
    gr->SetLineWidth(2);
    gr->SetLineColor(color);
    gr->SetMarkerColor(color);
    gr->Draw("P SAME");
    return gr;
} //enddef

} // anon

// ------------------------------------------------------------
// Main driver
// ------------------------------------------------------------
void compute_bin_centered_cross_sections(
    const std::vector<Binning>& binning_scheme,
    const std::string& radcorr_xsec_json_dir,
    const std::string& output_dir,
    int n_steps,
    const ModelPaths& paths,
    bool vgg_globalfit,
    ModelChoice model_choice
) {
    std::cout << "[bincenter] Starting bin-centering corrections..." << std::endl;
    std::cout << "[bincenter] Model choice: ";
    switch (model_choice) {
        case ModelChoice::Both: std::cout << "Both models (averaged)"; break;
        case ModelChoice::VGGOnly: std::cout << "VGG only"; break;
        case ModelChoice::KM15Only: std::cout << "KM15 only"; break;
    }
    std::cout << std::endl;
    
    std::cout << "[bincenter] Creating output directories..." << std::endl;
    
    fs::create_directories(output_dir);
    fs::create_directories(fs::path(output_dir)/"jsons");
    fs::create_directories(fs::path(output_dir)/"plots");

    const auto xB_bins = uniqueRanges(binning_scheme, 'x');
    const auto Q2_all  = uniqueRanges(binning_scheme, 'Q');
    const auto t_all   = uniqueRanges(binning_scheme, 't');
    const auto PHI_DEG = phiCentersDeg();

    // Keep 10.59 in the list; we will silently skip if it is not available yet.
    const std::vector<std::string> energies = {"10.59","10.60","10.2"};

    int total_energy_count = energies.size();
    int current_energy_index = 0;

    for (const auto& E : energies) {
        current_energy_index++;
        std::cout << "[bincenter] ===========================================" << std::endl;
        std::cout << "[bincenter] Processing energy E = " << E << " (" 
                  << current_energy_index << "/" << total_energy_count << ")" << std::endl;

        const std::string rc_path =
            (fs::path(radcorr_xsec_json_dir) / ("rad_corrected_xsec_" + E + ".json")).string();

        // Quietly skip if file is missing
        if (!fs::exists(rc_path)) {
            std::cout << "[bincenter] Skipping E=" << E << " (no RC file yet)\n";
            continue;
        } //endif

        // Load and validate
        std::cout << "[bincenter] Loading radiative corrections from: " << rc_path << std::endl;
        json j_rc = load_json(rc_path);
        if (j_rc.is_null() || j_rc.empty() || !j_rc.contains("bins")) {
            std::cout << "[bincenter] Skipping E=" << E << " (RC file present but malformed)\n";
            continue;
        } //endif

        // Output accumulator for this energy
        json jout;
        jout["energy"] = E;
        jout["bins"]   = json::object();

        int total_xb_bins = xB_bins.size();
        int current_xb_index = 0;

        // Loop the physics bins
        for (int ix = 0; ix < total_xb_bins; ++ix) {
            current_xb_index++;
            const auto xb = xB_bins[ix];

            std::cout << "[bincenter]   Processing xB bin " << current_xb_index << "/" << total_xb_bins 
                      << ": [" << xb.first << ", " << xb.second << "]" << std::endl;

            // Slices that actually occur for this xB in the scheme
            std::set<std::pair<double,double>> qs, ts;
            for (const auto& b : binning_scheme) {
                if (std::make_pair(b.xBmin,b.xBmax) == xb) {
                    qs.emplace(b.Q2min,b.Q2max);
                    ts.emplace(b.tmin,b.tmax);
                } //endif
            } //endfor
            std::vector<std::pair<double,double>> Q2_slice(qs.begin(), qs.end());
            std::vector<std::pair<double,double>> t_slice(ts.begin(),  ts.end());

            int total_q2_bins = Q2_slice.size();
            int total_t_bins = t_slice.size();
            int current_bin_count = 0;
            int total_bins_for_xb = total_q2_bins * total_t_bins;

            std::cout << "[bincenter]     Found " << total_q2_bins << " Q2 bins and " 
                      << total_t_bins << " t bins (" << total_bins_for_xb << " total bins)" << std::endl;

            for (const auto& q2r : Q2_slice) {
                const int iQ_global = findIndex(q2r, Q2_all);
                if (iQ_global < 0) continue;
                for (const auto& tr : t_slice) {
                    current_bin_count++;
                    const int it_global = findIndex(tr, t_all);
                    if (it_global < 0) continue;

                    const std::string bkey =
                        "(" + std::to_string(ix) + "," + std::to_string(iQ_global) + "," + std::to_string(it_global) + ")";

                    if (!j_rc["bins"].contains(bkey)) {
                        continue;
                    }
                    const json& cell_in = j_rc["bins"][bkey];
                    if (!has_bc_cell(cell_in)) {
                        continue;
                    }

                    auto do_one_hel = [&](const char* node, Helicity hel)->json{
                        json out;
                        if (!cell_in.contains(node)) {
                            return out;
                        }
                        const json& h = cell_in[node];
                        if (!has_bins12(h)) {
                            return out;
                        }

                        const auto& phi = h["phi"];
                        const auto& xs  = h["xsec"];
                        const auto& xe  = h["xsec_err"];

                        std::vector<double> phi_bc(N_PHI_BINS), y(N_PHI_BINS), e(N_PHI_BINS);
                        std::vector<double> fbin_used(N_PHI_BINS, 1.0), fbin_km15(N_PHI_BINS, 1.0), fbin_vgg(N_PHI_BINS, 1.0), fbin_sys(N_PHI_BINS, 0.0);

                        for (int ip=0; ip<N_PHI_BINS; ++ip){
                            const double phi_c = jgetd(phi, ip, PHI_DEG[ip]);
                            const auto   phi_edges = phiEdgesDegForBin(ip);

                            double fkm = 0.0, fvg = 0.0;
                            compute_fbin_for_point(
                                xb.first, xb.second,
                                q2r.first, q2r.second,
                                tr.first,  tr.second,
                                phi_c, phi_edges,
                                n_steps,
                                std::stod(E), hel,
                                paths, vgg_globalfit,
                                model_choice,
                                fkm, fvg
                            );

                            // Choose final factor based on model choice
                            double f_final = 1.0;
                            double f_sys   = 0.0;
                            
                            switch (model_choice) {
                                case ModelChoice::Both:
                                    if (finite_pos(fkm) && finite_pos(fvg)) {
                                        f_final = 0.5*(fkm + fvg);
                                        f_sys   = std::fabs(fkm - fvg);
                                    } else if (finite_pos(fkm)) {
                                        f_final = fkm;
                                        f_sys   = 1.0;
                                    } else if (finite_pos(fvg)) {
                                        f_final = fvg;
                                        f_sys   = 1.0;
                                    }
                                    break;
                                case ModelChoice::VGGOnly:
                                    if (finite_pos(fvg)) {
                                        f_final = fvg;
                                        f_sys   = 0.0; // No systematic when using single model
                                    }
                                    break;
                                case ModelChoice::KM15Only:
                                    if (finite_pos(fkm)) {
                                        f_final = fkm;
                                        f_sys   = 0.0; // No systematic when using single model
                                    }
                                    break;
                            }

                            const double x  = jgetd(xs, ip, 0.0);
                            const double ex = std::max(0.0, jgetd(xe, ip, 0.0));

                            phi_bc[ip]   = phi_c;
                            y[ip]        = f_final * x;
                            e[ip]        = f_final * ex;  // scale statistical error by same factor
                            fbin_used[ip]= f_final;
                            fbin_km15[ip]= (finite_pos(fkm) ? fkm : 0.0);
                            fbin_vgg[ip] = (finite_pos(fvg) ? fvg : 0.0);
                            fbin_sys[ip] = f_sys;
                        } //endfor

                        out["phi"]       = phi_bc;
                        out["xsec"]      = y;
                        out["xsec_err"]  = e;
                        out["fbin_used"] = fbin_used;
                        out["fbin_km15"] = fbin_km15;
                        out["fbin_vgg"]  = fbin_vgg;
                        out["fbin_sys"]  = fbin_sys;
                        return out;
                    }; //enddef

                    json hp_bc = do_one_hel("helicity_plus",  Helicity::Plus);
                    json hm_bc = do_one_hel("helicity_minus", Helicity::Minus);
                    if (hp_bc.is_null() && hm_bc.is_null()) {
                        continue;
                    }

                    jout["bins"][bkey] = {
                        {"helicity_plus",  hp_bc},
                        {"helicity_minus", hm_bc}
                    };
                } //endfor
            } //endfor
        } //endfor

        // If nothing was produced for this energy, do not write files or make plots
        if (jout["bins"].empty()) {
            std::cout << "[bincenter] No valid bins at E=" << E << ", skipping save/plots\n";
            continue;
        } //endif

        // Save JSON for this energy
        {
            const std::string out_json =
                (fs::path(output_dir)/"jsons"/("bin_centered_xsec_"+E+".json")).string();
            std::cout << "[bincenter] Saving results to: " << out_json << std::endl;
            save_json(jout, out_json);
            std::cout << "[bincenter] Successfully wrote " << out_json << "\n";
        }

        // Plot helper
        auto makeCanvas = [&](bool overlay){
            std::cout << "[bincenter] Generating " << (overlay ? "overlay" : "corrected-only") << " plots for E=" << E << std::endl;
            for (int ix = 0; ix < (int)xB_bins.size(); ++ix) {
                const auto xb = xB_bins[ix];

                std::set<std::pair<double,double>> qs, ts;
                for (const auto& b : binning_scheme) {
                    if (std::make_pair(b.xBmin,b.xBmax) == xb) {
                        qs.emplace(b.Q2min,b.Q2max);
                        ts.emplace(b.tmin,b.tmax);
                    } //endif
                } //endfor
                std::vector<std::pair<double,double>> Q2_slice(qs.begin(), qs.end());
                std::vector<std::pair<double,double>> t_slice(ts.begin(), ts.end());
                if (Q2_slice.empty() || t_slice.empty()) continue;

                const int nrows = (int)t_slice.size();
                const int ncols = (int)Q2_slice.size();
                const int W = 280*ncols + 160;
                const int H = 240*nrows + 170;

                std::ostringstream cname;
                cname << "c_bincenter_"<<E<<"_xB"<<ix<<(overlay ? "_overlay" : "_bc");
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
                head.SetTextSize(0.36);
                std::ostringstream tit;
                tit << (overlay ? "Before vs. After Bin-Centering d#sigma/d#phi"
                                : "Bin-Centered d#sigma/d#phi");
                tit << Form("  %s GeV   x_{B} in [%.2g, %.2g]", E.c_str(), xb.first, xb.second);
                head.DrawLatex(0.5, 0.55, tit.str().c_str());

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

                        TH1* frame = gPad->DrawFrame(0.0, 1e-4, 360.0, 1.0);
                        frame->GetXaxis()->SetLabelSize(0.0001);
                        frame->GetXaxis()->SetTitle("#phi (deg)");
                        frame->GetYaxis()->SetTitle(overlay ? "d#sigma/d#phi" : "d#sigma/d#phi (bin-centered)");
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

                        if (!jout["bins"].contains(bkey)) continue;
                        const json& jb_bc = jout["bins"][bkey];

                        const bool have_before = overlay && j_rc["bins"].contains(bkey);
                        const json& jb_before = have_before ? j_rc["bins"][bkey] : json::object();

                        // draw order: before first, after second
                        TGraphErrors* gup = nullptr;
                        TGraphErrors* gum = nullptr;
                        if (have_before) {
                            gup = draw_hel_curve(jb_before["helicity_plus"],  24, kGray+2, frame);
                            gum = draw_hel_curve(jb_before["helicity_minus"], 26, kGray+2, frame);
                        } //endif
                        TGraphErrors* gcp = draw_hel_curve(jb_bc["helicity_plus"],  20, kBlue+1, frame);
                        TGraphErrors* gcm = draw_hel_curve(jb_bc["helicity_minus"], 25, kRed+1,  frame);
                    } //endfor
                } //endfor

                const std::string outP =
                    (fs::path(output_dir)/"plots"/(std::string("bin_centered_xsec_")+E+"_xB_"+std::to_string(ix)+(overlay ? "_overlay.png" : ".png"))).string();
                c->SaveAs(outP.c_str());
                delete c;
            } //endfor
        }; //enddef

        makeCanvas(false); // after BC only
        makeCanvas(true);  // before vs after BC
        
        std::cout << "[bincenter] Completed processing for E = " << E << std::endl;
    } //endfor

    std::cout << "[bincenter] ===========================================" << std::endl;
    std::cout << "[bincenter] Finished bin-centering correction." << std::endl;
} //enddef