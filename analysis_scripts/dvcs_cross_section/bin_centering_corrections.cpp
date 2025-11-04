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

// OpenMP for parallel processing (hard-capped to 5 threads)
#ifdef _OPENMP
#include <omp.h>
#endif

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

// ---------- small helpers ----------
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
static inline bool   finite_pos(double v){ return std::isfinite(v) && v>0.0; } //enddef

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
    double vgg_center  = 0.0;

    if (model_choice == ModelChoice::Both || model_choice == ModelChoice::KM15Only) {
        km15_center = km15_xs(xB_c, Q2_c, tpos_c, phi_center_deg, Ebeam, hel, paths);
    }
    if (model_choice == ModelChoice::Both || model_choice == ModelChoice::VGGOnly) {
        vgg_center  = vgg_xs(xB_c, Q2_c, tpos_c, phi_center_deg, Ebeam, hel, paths, vgg_globalfit);
    }

    const bool need_km15 = (model_choice == ModelChoice::Both || model_choice == ModelChoice::KM15Only) && finite_pos(km15_center);
    const bool need_vgg  = (model_choice == ModelChoice::Both || model_choice == ModelChoice::VGGOnly)   && finite_pos(vgg_center);
    if (!need_km15 && !need_vgg) return;

    // sub-binning grids (uniform)
    const int nx = std::max(2, n_steps);
    const int nQ = std::max(2, n_steps);
    const int nt = std::max(2, n_steps);
    const int np = std::max(2, n_steps);

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

    double sum_km15 = 0.0; int cnt_km15 = 0;
    double sum_vgg  = 0.0; int cnt_vgg  = 0;

    for (double xb : xs){
        for (double q2 : Qs){
            for (double tp : ts){
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

    if (need_km15 && cnt_km15 > 0) {
        const double avg_km15 = sum_km15 / double(cnt_km15);
        if (finite_pos(avg_km15)) fbin_km15 = km15_center / avg_km15;
    }
    if (need_vgg && cnt_vgg > 0) {
        const double avg_vgg = sum_vgg / double(cnt_vgg);
        if (finite_pos(avg_vgg)) fbin_vgg = vgg_center / avg_vgg;
    }
} //enddef

// ---------------- plotting helpers (shared with plotter) ----------------
static TGraphErrors* draw_hel_curve(const json& h, int mstyle, int color){
    if (!has_bins12(h)) return nullptr; //endif
    const auto& xp = h["phi"];
    const auto& yp = h["xsec"];
    const auto& ep = h["xsec_err"];
    std::vector<double> x(N_PHI_BINS), y(N_PHI_BINS), e(N_PHI_BINS);
    for (int i=0;i<N_PHI_BINS;++i){
        x[i]=xp[i].get<double>();
        y[i]=yp[i].get<double>();
        e[i]=std::max(1e-12, ep[i].get<double>());
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

// For log-scale y-range across a whole canvas (like rad_corrected_cross_section.cpp)
static std::pair<double, double> calculateYRangeForCanvas(const json& jb_after_bins,
                                                         const json& jb_before_bins,
                                                         const std::vector<std::string>& bin_keys,
                                                         bool overlay) {
    double global_min = 1e10;
    double global_max = -1e10;

    auto scanHel = [&](const json& H){
        if (!H.contains("xsec") || !H.contains("xsec_err")) return;
        const auto& yp = H["xsec"];
        const auto& ep = H["xsec_err"];
        for (int i=0;i<N_PHI_BINS;++i){
            const double y = jgetd(yp, i, 0.0);
            const double e = jgetd(ep, i, 0.0);
            if (y > 0.0) {
                global_min = std::min(global_min, std::max(1e-10, y - e));
                global_max = std::max(global_max, y + e);
            }
        }
    }; //enddef

    for (const auto& bkey : bin_keys) {
        if (jb_after_bins.contains(bkey)) {
            const json& after = jb_after_bins[bkey];
            if (after.contains("helicity_plus"))  scanHel(after["helicity_plus"]);
            if (after.contains("helicity_minus")) scanHel(after["helicity_minus"]);
        }
        if (overlay && jb_before_bins.contains(bkey)) {
            const json& before = jb_before_bins[bkey];
            if (before.contains("helicity_plus"))  scanHel(before["helicity_plus"]);
            if (before.contains("helicity_minus")) scanHel(before["helicity_minus"]);
        }
    } //endfor

    if (global_max <= 0) {
        global_min = 1e-4;
        global_max = 1.0;
    } else {
        double ymin10 = std::pow(10.0, std::floor(std::log10(global_min)));
        double ymax10 = std::pow(10.0, std::ceil (std::log10(global_max)));
        global_min = std::max(1e-4, ymin10*0.5);
        global_max = ymax10*2.0;
    }
    return {global_min, global_max};
} //enddef

} // anon

// ============================================================
// COMPUTE-ONLY: writes bin_centered_xsec_<E>.json (no plots)
// ============================================================
void compute_bin_centered_cross_sections(
    const std::vector<Binning>& binning_scheme,
    const std::string& radcorr_xsec_json_dir,  // where rad_corrected_xsec_<E>.json lives (input)
    const std::string& output_dir,             // base output dir; JSONs go to output_dir/jsons
    int n_steps,
    const ModelPaths& paths,
    bool vgg_globalfit,
    ModelChoice model_choice
) {
    std::cout << "[bincenter] Starting bin-centering computation..." << std::endl;
    std::cout << "[bincenter] Model choice = "
              << (model_choice==ModelChoice::Both? "Both" : (model_choice==ModelChoice::VGGOnly? "VGGOnly" : "KM15Only"))
              << ", n_steps = " << n_steps << std::endl;

    #ifdef _OPENMP
    int hard_cap = 5;
    int want     = omp_get_max_threads();
    int use      = std::min(hard_cap, std::max(1, want));
    omp_set_num_threads(use);
    std::cout << "[bincenter] OpenMP enabled with " << use << " worker(s) (hard-capped at 5)\n";
    #else
    std::cout << "[bincenter] OpenMP not available; running single-threaded\n";
    #endif

    fs::create_directories(output_dir);
    fs::create_directories(fs::path(output_dir)/"jsons");

    const auto xB_bins = uniqueRanges(binning_scheme, 'x');
    const auto Q2_all  = uniqueRanges(binning_scheme, 'Q');
    const auto t_all   = uniqueRanges(binning_scheme, 't');
    const auto PHI_DEG = phiCentersDeg();

    // energies we attempt (skip missing quietly)
    const std::vector<std::string> energies = {"10.59","10.60","10.2"};

    int e_idx = 0;
    for (const auto& E : energies) {
        ++e_idx;
        const std::string rc_path = (fs::path(radcorr_xsec_json_dir) / ("rad_corrected_xsec_" + E + ".json")).string();
        if (!fs::exists(rc_path)) {
            std::cout << "[bincenter] Skip E=" << E << " (missing " << rc_path << ")\n";
            continue;
        }
        std::cout << "[bincenter] ===========================================\n";
        std::cout << "[bincenter] Energy " << E << " (" << e_idx << "/" << energies.size() << ") -> " << rc_path << "\n";

        json j_rc = load_json(rc_path);
        if (j_rc.is_null() || j_rc.empty() || !j_rc.contains("bins")) {
            std::cout << "[bincenter] JSON malformed for E=" << E << " (no bins). Skipping.\n";
            continue;
        }

        json jout; // result for this energy
        jout["energy"] = E;
        jout["bins"]   = json::object();

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
            std::vector<std::pair<double,double>> t_slice(ts.begin(),  ts.end());

            // Gather (Q2,t) tasks
            std::vector<std::pair<std::pair<double,double>, std::pair<double,double>>> tasks;
            tasks.reserve(Q2_slice.size()*t_slice.size());
            for (const auto& q2r : Q2_slice) {
                for (const auto& tr : t_slice) {
                    tasks.emplace_back(q2r, tr);
                }
            }

            #pragma omp parallel for schedule(dynamic)
            for (size_t it = 0; it < tasks.size(); ++it) {
                const auto q2r = tasks[it].first;
                const auto tr  = tasks[it].second;

                const int iQ_global = findIndex(q2r, Q2_all);
                const int it_global = findIndex(tr,  t_all);
                if (iQ_global < 0 || it_global < 0) continue;

                const std::string bkey =
                    "(" + std::to_string(ix) + "," + std::to_string(iQ_global) + "," + std::to_string(it_global) + ")";

                if (!j_rc["bins"].contains(bkey)) continue;
                const json& cell_in = j_rc["bins"][bkey];
                if (!has_bc_cell(cell_in)) continue;

                auto do_one_hel = [&](const char* node, Helicity hel)->json{
                    json out;
                    if (!cell_in.contains(node)) return out;
                    const json& h = cell_in[node];
                    if (!has_bins12(h)) return out;

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

                        // determine final factor
                        double f_final = 1.0;
                        double f_sys   = 0.0;
                        switch (model_choice) {
                            case ModelChoice::Both:
                                if (finite_pos(fkm) && finite_pos(fvg)) {
                                    f_final = 0.5*(fkm + fvg);
                                    f_sys   = std::fabs(fkm - fvg);
                                } else if (finite_pos(fkm)) {
                                    f_final = fkm; f_sys = 1.0;
                                } else if (finite_pos(fvg)) {
                                    f_final = fvg; f_sys = 1.0;
                                }
                                break;
                            case ModelChoice::VGGOnly:
                                if (finite_pos(fvg)) f_final = fvg;
                                break;
                            case ModelChoice::KM15Only:
                                if (finite_pos(fkm)) f_final = fkm;
                                break;
                        } //endswitch

                        const double x  = jgetd(xs, ip, 0.0);
                        const double ex = std::max(0.0, jgetd(xe, ip, 0.0));

                        phi_bc[ip]   = phi_c;
                        y[ip]        = f_final * x;
                        e[ip]        = f_final * ex;  // scale stat error
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
                if (hp_bc.is_null() && hm_bc.is_null()) continue;

                #ifdef _OPENMP
                #pragma omp critical
                #endif
                {
                    jout["bins"][bkey] = {
                        {"helicity_plus",  hp_bc},
                        {"helicity_minus", hm_bc}
                    };
                }
            } //end parallel for
        } //endfor xB

        // Write out JSON for this energy
        const std::string out_json = (fs::path(output_dir)/"jsons"/("bin_centered_xsec_"+E+".json")).string();
        save_json(jout, out_json);
        std::cout << "[bincenter] Wrote " << out_json << "\n";
    } //endfor energies

    std::cout << "[bincenter] Finished bin-centering computation.\n";
} //enddef

// ============================================================
// PLOT-ONLY: reads bin_centered_xsec_<E>.json (+ before BC)
// makes corrected-only and overlay plots
// ============================================================
void plot_bin_centered_cross_sections(
    const std::vector<Binning>& binning_scheme,
    const std::string& radcorr_xsec_json_dir, // input: before-BC "rad_corrected_xsec_<E>.json"
    const std::string& bincenter_json_dir,    // input: after-BC  "bin_centered_xsec_<E>.json"
    const std::string& plots_output_dir       // output: plots go here
) {
    fs::create_directories(plots_output_dir);

    const auto xB_bins = uniqueRanges(binning_scheme, 'x');
    const auto Q2_all  = uniqueRanges(binning_scheme, 'Q');
    const auto t_all   = uniqueRanges(binning_scheme, 't');

    const std::vector<std::string> energies = {"10.59","10.60","10.2"};

    auto makeCanvasSet = [&](const std::string& E, const json& j_before, const json& j_after, bool overlay){
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
            std::vector<std::pair<double,double>> t_slice(ts.begin(),  ts.end());
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

            // Title (match rad_corrected_cross_section.cpp sizing/position)
            pTop->cd();
            TLatex head;
            head.SetNDC(); head.SetTextAlign(22);
            head.SetTextFont(42);
            head.SetTextSize(0.22);
            std::ostringstream tit;
            tit << (overlay ? "unc vs. bin-centered d#sigma/d#phi"
                            : "bin-centered d#sigma/d#phi");
            tit << Form("  %s GeV x_{B} #in [%.2g, %.2g]", E.c_str(), xb.first, xb.second);
            head.DrawLatex(0.5, 0.55, tit.str().c_str());

            // First pass: collect all bin keys for this canvas to determine y-range
            std::vector<std::string> canvas_bin_keys;
            for (int r=0; r<nrows; ++r) {
                const int it_global = findIndex(t_slice[r], t_all);
                if (it_global < 0) continue;
                for (int cc=0; cc<ncols; ++cc) {
                    const int iQ_global = findIndex(Q2_slice[cc], Q2_all);
                    if (iQ_global < 0) continue;
                    const std::string bkey =
                        "(" + std::to_string(ix) + "," + std::to_string(iQ_global) + "," + std::to_string(it_global) + ")";
                    canvas_bin_keys.push_back(bkey);
                }
            }

            // Compute common y-range for the canvas
            auto [ymin_canvas, ymax_canvas] = calculateYRangeForCanvas(
                j_after["bins"],
                j_before.contains("bins") ? j_before["bins"] : json::object(),
                canvas_bin_keys,
                overlay
            );

            // Draw panels
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

                    // after-BC (required)
                    if (!j_after["bins"].contains(bkey)) continue;
                    const json& jb_after = j_after["bins"][bkey];

                    // before-BC (optional overlay)
                    const bool   have_before = overlay && j_before.contains("bins") && j_before["bins"].contains(bkey);
                    const json&  jb_before  = (have_before ? j_before["bins"][bkey] : json::object());

                    // draw order: before first (grey open), after second (colored filled)
                    TGraphErrors* gup = nullptr;
                    TGraphErrors* gum = nullptr;
                    if (have_before) {
                        if (jb_before.contains("helicity_plus"))
                            gup = draw_hel_curve(jb_before["helicity_plus"],  24 /*open circle*/,   kGray+2);
                        if (jb_before.contains("helicity_minus"))
                            gum = draw_hel_curve(jb_before["helicity_minus"], 26 /*open triangle*/, kGray+2);
                    }
                    TGraphErrors* gcp = nullptr;
                    TGraphErrors* gcm = nullptr;
                    if (jb_after.contains("helicity_plus"))
                        gcp = draw_hel_curve(jb_after["helicity_plus"],  20 /*filled circle*/,  kBlue+1);
                    if (jb_after.contains("helicity_minus"))
                        gcm = draw_hel_curve(jb_after["helicity_minus"], 25 /*filled square*/,  kRed+1);

                    // per-panel label (Q2,t) — match sizes/placement from rad code
                    TLatex lab;
                    lab.SetNDC(); lab.SetTextSize(0.07); lab.SetTextAlign(11);
                    lab.SetTextFont(42);
                    lab.DrawLatex(0.15, 0.94,
                        Form("Q^{2} #in [%.2g, %.2g], -t #in [%.2g, %.2g]",
                             Q2_slice[cc].first, Q2_slice[cc].second,
                             t_slice[r].first,   t_slice[r].second));

                    // bottom-left legend (match rad code)
                    if (gcp || gcm || gup || gum){
                        double x1=0.16, y1=0.18, x2=0.68, y2=0.42;
                        TLegend* leg = new TLegend(x1,y1,x2,y2);
                        leg->SetBorderSize(1);
                        leg->SetLineColor(kBlack);
                        leg->SetFillColor(kWhite);
                        leg->SetFillStyle(1001);
                        leg->SetTextFont(42);
                        leg->SetTextSize(0.048);
                        if (gcp) leg->AddEntry(gcp, "bin-centered  + helicity", "lep");
                        if (gcm) leg->AddEntry(gcm, "bin-centered  - helicity", "lep");
                        if (gup) leg->AddEntry(gup, "before BC     + helicity", "lep");
                        if (gum) leg->AddEntry(gum, "before BC     - helicity", "lep");
                        leg->Draw();
                    } //endif
                } //endfor cols
            } //endfor rows

            const std::string outP =
                (fs::path(plots_output_dir) / (std::string("bin_centered_xsec_") + E + "_xB_" + std::to_string(ix) + (overlay ? "_overlay.png" : ".png"))).string();
            c->SaveAs(outP.c_str());
            delete c;
        } //endfor xB
    }; //end lambda makeCanvasSet

    for (const auto& E : energies) {
        const std::string before_path = (fs::path(radcorr_xsec_json_dir) / ("rad_corrected_xsec_" + E + ".json")).string();
        const std::string after_path  = (fs::path(bincenter_json_dir)    / ("bin_centered_xsec_" + E + ".json")).string();

        if (!fs::exists(after_path)) {
            std::cout << "[bincenter-plot] Skip E=" << E << " (missing " << after_path << ")\n";
            continue;
        }

        json j_before = load_json(before_path); // may be empty if missing; overlay will just omit
        json j_after  = load_json(after_path);

        if (j_after.is_null() || !j_after.contains("bins")) {
            std::cout << "[bincenter-plot] Malformed after-BC JSON for E=" << E << "\n";
            continue;
        }

        std::cout << "[bincenter-plot] Plotting E=" << E << "\n";
        makeCanvasSet(E, j_before, j_after, false); // corrected-only
        makeCanvasSet(E, j_before, j_after, true);  // before vs after
    } //endfor energies

    std::cout << "[bincenter-plot] Finished plotting bin-centered cross sections.\n";
} //enddef