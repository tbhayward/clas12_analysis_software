// models_vs_data_plots.cpp
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
#include <atomic>
#include <chrono>
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

#ifdef _OPENMP
#include <omp.h>
#endif

#include <nlohmann/json.hpp>
using json = nlohmann::json;
namespace fs = std::filesystem;

namespace {

constexpr int N_PHI_BINS = 12;

// ------------ style bootstrap ------------
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

// Force line-buffered std::cout so progress prints immediately.
struct CoutUnitbuf {
    CoutUnitbuf() { std::cout.setf(std::ios::unitbuf); }
} _cout_unitbuf_guard;

// ---------- small helpers ----------
static inline void log_info(const std::string& msg) {
    std::cout << "[models] " << msg << std::endl;
}
static inline void log_warn(const std::string& msg) {
    std::cout << "[models][warn] " << msg << std::endl;
}

#ifdef _OPENMP
static inline int configure_omp_workers(int hard_cap = 5) {
    // eliminate surprises that cause slowdowns / nested teams
    omp_set_dynamic(0);
    omp_set_nested(0);
    omp_set_max_active_levels(1);
    int sys_cap = omp_get_max_threads();
    int use = std::min(hard_cap, std::max(1, sys_cap));
    omp_set_num_threads(use);
    log_info("OpenMP enabled with " + std::to_string(use) + " worker(s) (hard-capped at 5).");
    return use;
}
#else
static inline int configure_omp_workers(int /*hard_cap*/ = 5) {
    log_info("OpenMP not available; running single-threaded.");
    return 1;
}
#endif

static inline std::vector<double> phiCentersDeg() {
    std::vector<double> d(N_PHI_BINS);
    const double step = 360.0 / double(N_PHI_BINS);
    for (int i=0; i<N_PHI_BINS; ++i) d[i] = (i+0.5)*step;
    return d;
}

static void drawDegreeTicks(double xmin, double ymin, double xmax, double labelSize) {
    TGaxis* ax = new TGaxis(xmin, ymin, xmax, ymin, 0.0, 360.0, 4, "");
    ax->SetLabelFont(42);
    ax->SetLabelSize(labelSize);
    ax->SetLabelOffset(0.012);
    ax->SetTitle("");
    ax->SetTickSize(0.02);
    ax->Draw();
}

static std::vector<std::pair<double,double>> uniqueRanges(const std::vector<Binning>& scheme, char which) {
    std::set<std::pair<double,double>> s;
    for (const auto& b : scheme) {
        if      (which == 'x') s.emplace(b.xBmin, b.xBmax);
        else if (which == 'Q') s.emplace(b.Q2min, b.Q2max);
        else if (which == 't') s.emplace(b.tmin, b.tmax);
    }
    return std::vector<std::pair<double,double>>(s.begin(), s.end());
}

static inline int findIndex(const std::pair<double,double>& range,
                            const std::vector<std::pair<double,double>>& ranges) {
    for (int i=0; i<(int)ranges.size(); ++i) if (ranges[i] == range) return i;
    return -1;
}

static json load_json(const std::string& path) {
    std::ifstream f(path);
    if (!f.is_open()) {
        log_warn("Failed to open " + path);
        return json();
    }
    json j; f >> j; return j;
}
static void save_json(const json& j, const std::string& path) {
    fs::create_directories(fs::path(path).parent_path());
    std::ofstream ofs(path);
    if (!ofs.is_open()) {
        log_warn("Failed to write " + path);
        return;
    }
    ofs << std::setw(2) << j << "\n";
}

static inline bool has_bins12(const json& h) {
    return h.contains("phi") && h.contains("xsec") && h.contains("xsec_err")
        && h["phi"].size()==N_PHI_BINS && h["xsec"].size()==N_PHI_BINS && h["xsec_err"].size()==N_PHI_BINS;
}
static inline bool positive(double v) { return std::isfinite(v) && v > 0.0; }
static inline double jgetd(const json& a, size_t i, double def=0.0) {
    try { return a[i].get<double>(); } catch(...) { return def; }
}
static inline double midpoint(double a, double b) { return 0.5*(a+b); }

// Draw the data points (bin-centered cross sections)
static TGraphErrors* draw_data_hel(const json& h, int mstyle, int color, TH1* /*frame*/, double& ymax_accum) {
    if (!has_bins12(h)) return nullptr;
    const auto& xp = h["phi"];
    const auto& yp = h["xsec"];
    const auto& ep = h["xsec_err"];
    std::vector<double> x(N_PHI_BINS), y(N_PHI_BINS), e(N_PHI_BINS);
    for (int i=0; i<N_PHI_BINS; ++i) {
        x[i] = xp[i].get<double>();
        y[i] = yp[i].get<double>();
        e[i] = std::max(1e-12, ep[i].get<double>());
        ymax_accum = std::max(ymax_accum, y[i] + e[i]);
    }
    auto* gr = new TGraphErrors(N_PHI_BINS, x.data(), y.data(), nullptr, e.data());
    gr->SetMarkerStyle(mstyle);
    gr->SetMarkerSize(1.0);
    gr->SetLineWidth(2);
    gr->SetLineColor(color);
    gr->SetMarkerColor(color);
    gr->Draw("P SAME");
    return gr;
}

// Envelope helper for y-range from predictions in one panel
static void scan_prediction_envelope(const json& pb, double& mmin, double& mmax) {
    auto upd = [&](const char* key) {
        if (!pb.contains(key)) return;
        const auto& v = pb[key];
        for (size_t i=0; i<v.size(); ++i) {
            double val = 0.0; try { val = v[i].get<double>(); } catch(...) { val = 0.0; }
            if (positive(val)) { mmin = std::min(mmin, val); mmax = std::max(mmax, val); }
        }
    };
    upd("bh");
    upd("vgg_plus");  upd("vgg_minus");
    upd("km15_plus"); upd("km15_minus");
}

// Compute canvas-wide y-range including data and prediction envelopes
static std::pair<double,double> computeYRangeCanvas(
    const json& data_bins,               // bin-centered data
    const json& pred_bins,               // model predictions
    const std::vector<std::string>& bks  // bin keys in this canvas
) {
    double global_min = 1e10;
    double global_max = -1e10;

    auto collectData = [&](const json& h) {
        if (!h.contains("xsec") || !h.contains("xsec_err")) return;
        const auto& yp = h["xsec"];
        const auto& ep = h["xsec_err"];
        for (int i=0; i<N_PHI_BINS; ++i) {
            double y = jgetd(yp, i, 0.0);
            double e = jgetd(ep, i, 0.0);
            if (y > 0.0) {
                global_min = std::min(global_min, std::max(1e-12, y - e));
                global_max = std::max(global_max, y + e);
            }
        }
    };

    for (const auto& bkey : bks) {
        if (data_bins.contains(bkey)) {
            const json& jb = data_bins[bkey];
            if (jb.contains("helicity_plus"))  collectData(jb["helicity_plus"]);
            if (jb.contains("helicity_minus")) collectData(jb["helicity_minus"]);
        }
        if (pred_bins.contains(bkey)) {
            double mmin = 1e10, mmax = -1e10;
            scan_prediction_envelope(pred_bins[bkey], mmin, mmax);
            if (mmax > 0.0) {
                global_min = std::min(global_min, mmin);
                global_max = std::max(global_max, mmax);
            }
        }
    }

    if (!(global_max > 0.0)) return {1e-4, 1.0};
    double ymin = std::pow(10.0, std::floor(std::log10(std::max(1e-12, global_min))));
    ymin = std::max(1e-4, ymin*0.5);
    double ymax = std::pow(10.0, std::ceil(std::log10(global_max)));
    ymax *= 2.0;
    return {ymin, ymax};
}

} // anon

// ============================================================
// 1) HEAVY STEP: compute and save model predictions JSON
// ============================================================
void compute_model_predictions(
    const std::vector<Binning>& binning_scheme,
    const std::string& bincenter_json_dir,  // read bin_centered_xsec_<E>.json to know which bins exist
    const std::string& output_dir,          // write jsons/model_predictions_<E>.json here
    const ModelPaths& paths,
    bool vgg_globalfit,
    int phi_dense
) {
    log_info("Starting compute_model_predictions...");
    const int omp_workers = configure_omp_workers(5);
    if (phi_dense < 2) phi_dense = 2;

    fs::create_directories(output_dir);
    fs::create_directories(fs::path(output_dir)/"jsons");

    const auto xB_bins = uniqueRanges(binning_scheme, 'x');
    const auto Q2_all  = uniqueRanges(binning_scheme, 'Q');
    const auto t_all   = uniqueRanges(binning_scheme, 't');
    const auto PHI_DEG = phiCentersDeg();

    const std::vector<std::string> energies = {"10.59","10.60","10.2"};

    int eidx = 0;
    for (const auto& E : energies) {
        ++eidx;
        const std::string bc_path = (fs::path(bincenter_json_dir) / ("bin_centered_xsec_" + E + ".json")).string();
        if (!fs::exists(bc_path)) {
            log_warn("Missing bin-centered JSON for E=" + E + "; skipping predictions for this energy.");
            continue;
        }
        log_info("Energy " + E + " (" + std::to_string(eidx) + "/" + std::to_string(energies.size()) + ") loading bins from " + bc_path);
        json j_bc = load_json(bc_path);
        if (j_bc.empty() || !j_bc.contains("bins")) {
            log_warn("Malformed bin-centered JSON for E=" + E + " (no bins). Skipping.");
            continue;
        }

        json jout;
        jout["energy"]    = E;
        jout["phi_dense"] = phi_dense;
        jout["bins"]      = json::object();

        // Build panel task list for this energy
        struct Task {
            int ix;
            int iQ_global;
            int it_global;
            std::pair<double,double> xb;
            std::pair<double,double> q2r;
            std::pair<double,double> tr;
            std::string bkey;
        };
        std::vector<Task> tasks;

        for (int ix = 0; ix < (int)xB_bins.size(); ++ix) {
            const auto xb = xB_bins[ix];

            std::set<std::pair<double,double>> qs, ts;
            for (const auto& b : binning_scheme) {
                if (std::make_pair(b.xBmin,b.xBmax) == xb) {
                    qs.emplace(b.Q2min,b.Q2max);
                    ts.emplace(b.tmin,b.tmax);
                }
            }
            std::vector<std::pair<double,double>> Q2_slice(qs.begin(), qs.end());
            std::vector<std::pair<double,double>> t_slice(ts.begin(),  ts.end());

            for (const auto& q2r : Q2_slice) {
                for (const auto& tr : t_slice) {
                    int iQ_global = findIndex(q2r, Q2_all);
                    int it_global = findIndex(tr,  t_all);
                    if (iQ_global < 0 || it_global < 0) continue;
                    std::string bkey = "(" + std::to_string(ix) + "," + std::to_string(iQ_global) + "," + std::to_string(it_global) + ")";
                    if (!j_bc["bins"].contains(bkey)) continue; // skip panels that do not exist in data
                    tasks.push_back({ix, iQ_global, it_global, xb, q2r, tr, bkey});
                }
            }
        }

        const size_t total = tasks.size();
        log_info("E=" + E + " tasks: " + std::to_string(total) + " panels to compute.");
        if (total == 0) {
            log_warn("E=" + E + " has zero panels. Skipping.");
            continue;
        }

        std::atomic<size_t> done{0};
        const size_t tick = std::max<size_t>(1, total/10);
        auto t_energy_begin = std::chrono::steady_clock::now();

        // Parallelize across PANELS. Inside each panel, run phi loop sequentially to avoid nested OMP overhead.
        #ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic)
        #endif
        for (long it = 0; it < (long)total; ++it) {
            const Task task = tasks[it];
            const std::string& bkey = task.bkey;

            // START line for this panel (print once when a worker picks it up)
            #ifdef _OPENMP
            #pragma omp critical
            #endif
            {
                log_info("E=" + E + " START panel " + bkey +
                         " [" + std::to_string(it+1) + "/" + std::to_string(total) + "]");
            }

            const double xB_c   = midpoint(task.xb.first,   task.xb.second);
            const double Q2_c   = midpoint(task.q2r.first,  task.q2r.second);
            const double tpos_c = midpoint(task.tr.first,   task.tr.second);
            const double Ebeam  = std::stod(E);

            std::vector<double> ph(phi_dense);
            for (int i=0; i<phi_dense; ++i) {
                ph[i] = (phi_dense == 1) ? 0.0 : (360.0 * double(i) / double(phi_dense-1));
            }

            std::vector<double> y_bh(phi_dense, 0.0);
            std::vector<double> y_vg_p(phi_dense, 0.0), y_vg_m(phi_dense, 0.0);
            std::vector<double> y_km_p(phi_dense, 0.0), y_km_m(phi_dense, 0.0);

            // Sequential phi loop (fast enough; avoids nested parallel overhead)
            for (int i=0; i<phi_dense; ++i) {
                const double p = ph[i];
                // BH exact (helicity independent)
                double yb = 0.0;
                try { yb = vgg_bh_only(xB_c, Q2_c, tpos_c, p, Ebeam, paths, vgg_globalfit); } catch(...) { yb = 0.0; }
                y_bh[i] = positive(yb) ? yb : 0.0;

                // VGG
                double vg_p=0.0, vg_m=0.0;
                try {
                    vg_p = vgg_xs(xB_c, Q2_c, tpos_c, p, Ebeam, Helicity::Plus,  paths, vgg_globalfit);
                    vg_m = vgg_xs(xB_c, Q2_c, tpos_c, p, Ebeam, Helicity::Minus, paths, vgg_globalfit);
                } catch(...) { vg_p = vg_m = 0.0; }
                y_vg_p[i] = positive(vg_p) ? vg_p : 0.0;
                y_vg_m[i] = positive(vg_m) ? vg_m : 0.0;

                // KM15
                double km_p=0.0, km_m=0.0;
                try {
                    km_p = km15_xs(xB_c, Q2_c, tpos_c, p, Ebeam, Helicity::Plus,  paths);
                    km_m = km15_xs(xB_c, Q2_c, tpos_c, p, Ebeam, Helicity::Minus, paths);
                } catch(...) { km_p = km_m = 0.0; }
                y_km_p[i] = positive(km_p) ? km_p : 0.0;
                y_km_m[i] = positive(km_m) ? km_m : 0.0;
            }

            json pb;
            pb["meta"] = {
                {"xB",    {task.xb.first,  task.xb.second}},
                {"Q2",    {task.q2r.first, task.q2r.second}},
                {"tpos",  {task.tr.first,  task.tr.second}},
                {"centers", {{"xB", xB_c}, {"Q2", Q2_c}, {"tpos", tpos_c}}}
            };
            pb["phi_deg"]   = ph;
            pb["bh"]        = y_bh;
            pb["vgg_plus"]  = y_vg_p;
            pb["vgg_minus"] = y_vg_m;
            pb["km15_plus"] = y_km_p;
            pb["km15_minus"]= y_km_m;

            #ifdef _OPENMP
            #pragma omp critical
            #endif
            {
                jout["bins"][bkey] = pb;
            }

            size_t d = ++done;
            if (d % tick == 0 || d == total) {
                #ifdef _OPENMP
                #pragma omp critical
                #endif
                {
                    int pct = int(100.0 * double(d) / double(total));
                    log_info("E=" + E + " progress " + std::to_string(pct) +
                             "% (" + std::to_string(d) + "/" + std::to_string(total) + ")");
                }
            }
        } // panels loop

        const std::string out_json = (fs::path(output_dir)/"jsons"/("model_predictions_" + E + ".json")).string();
        save_json(jout, out_json);
        auto t_energy_end = std::chrono::steady_clock::now();
        auto ms_energy = std::chrono::duration_cast<std::chrono::milliseconds>(t_energy_end - t_energy_begin).count();
        log_info("E=" + E + " wrote " + out_json + " (" + std::to_string(ms_energy) + " ms).");
    } // energies

    log_info("Finished compute_model_predictions.");
}

// Add this helper near your other anon-namespace helpers (uses jgetd & scan_prediction_envelope already present)
static std::pair<double,double> computeYRangePanel(
    const json& data_bins,
    const json& pred_bins,
    const std::string& bkey
) {
    double global_min = 1e10;
    double global_max = -1e10;

    auto collectData = [&](const json& h) {
        if (!h.contains("xsec") || !h.contains("xsec_err")) return;
        const auto& yp = h["xsec"];
        const auto& ep = h["xsec_err"];
        for (int i=0; i<N_PHI_BINS; ++i) {
            double y = jgetd(yp, i, 0.0);
            double e = jgetd(ep, i, 0.0);
            if (y > 0.0) {
                global_min = std::min(global_min, std::max(1e-12, y - e));
                global_max = std::max(global_max, y + e);
            }
        }
    };

    if (data_bins.contains(bkey)) {
        const json& jb = data_bins[bkey];
        if (jb.contains("helicity_plus"))  collectData(jb["helicity_plus"]);
        if (jb.contains("helicity_minus")) collectData(jb["helicity_minus"]);
    }
    if (pred_bins.contains(bkey)) {
        double mmin = 1e10, mmax = -1e10;
        scan_prediction_envelope(pred_bins[bkey], mmin, mmax);
        if (mmax > 0.0) {
            global_min = std::min(global_min, mmin);
            global_max = std::max(global_max, mmax);
        }
    }

    if (!(global_max > 0.0)) return {1e-4, 1.0};

    // Log-scale friendly rounding with padding
    double ymin = std::pow(10.0, std::floor(std::log10(std::max(1e-12, global_min))));
    ymin = std::max(1e-4, ymin*0.5);
    double ymax = std::pow(10.0, std::ceil(std::log10(global_max)));
    ymax *= 2.0;
    return {ymin, ymax};
}



// ============================================================
// 2) FAST STEP: plot using saved predictions + bin-centered data
// ============================================================
static std::pair<double,double> panelYRangeLogSafe(const json& j_bc_cell,
                                                   const json* j_pred_cell) {
    auto scanHel = [](const json& H, double& gmin, double& gmax){
        if (!H.contains("xsec") || !H.contains("xsec_err")) return;
        const auto& y  = H["xsec"];
        const auto& ey = H["xsec_err"];
        for (int i=0;i<N_PHI_BINS;++i){
            const double v  = jgetd(y,  i, 0.0);
            const double dv = std::max(0.0, jgetd(ey, i, 0.0));
            if (v > 0.0) {
                gmin = std::min(gmin, std::max(1e-12, v - dv));
                gmax = std::max(gmax, v + dv);
            }
        }
    };

    double gmin = 1e10, gmax = -1e10;

    // data
    if (j_bc_cell.contains("helicity_plus"))  scanHel(j_bc_cell["helicity_plus"],  gmin, gmax);
    if (j_bc_cell.contains("helicity_minus")) scanHel(j_bc_cell["helicity_minus"], gmin, gmax);

    // predictions
    if (j_pred_cell){
        auto scanVec = [&](const char* key){
            if (!j_pred_cell->contains(key)) return;
            const auto& a = (*j_pred_cell)[key];
            for (size_t i=0;i<a.size();++i){
                const double v = jgetd(a, i, 0.0);
                if (positive(v)) { gmin = std::min(gmin, v); gmax = std::max(gmax, v); }
            }
        };
        scanVec("km15_plus");  scanVec("km15_minus");
        scanVec("vgg_plus");   scanVec("vgg_minus");
        scanVec("bh");
    }

    // fallbacks + log padding
    if (!(gmax > 0.0)) { gmin = 1e-4; gmax = 1.0; }
    double ymin10 = std::pow(10.0, std::floor(std::log10(std::max(1e-12, gmin))));
    double ymax10 = std::pow(10.0, std::ceil (std::log10(gmax)));
    gmin = std::max(1e-4, ymin10*0.8);
    gmax = ymax10*1.25;
    if (gmax <= gmin) gmax = gmin*10.0; // safety
    return {gmin, gmax};
}

void plot_models_vs_bincentered(
    const std::vector<Binning>& binning_scheme,
    const std::string& bincenter_json_dir,     // bin_centered_xsec_<E>.json (data)
    const std::string& predictions_json_dir,   // model_predictions_<E>.json (curves)
    const std::string& output_dir,             // PNGs under output_dir/plots
    int /*phi_dense_unused*/
) {
    log_info("Starting models vs. bin-centered plotting (using saved predictions)...");
    configure_omp_workers(5);

    fs::create_directories(output_dir);
    fs::create_directories(fs::path(output_dir)/"plots");

    const auto xB_bins = uniqueRanges(binning_scheme, 'x');
    const auto Q2_all  = uniqueRanges(binning_scheme, 'Q');
    const auto t_all   = uniqueRanges(binning_scheme, 't');

    log_info("Discovered unique bins: xB=" + std::to_string(xB_bins.size())
             + ", Q2=" + std::to_string(Q2_all.size())
             + ", t="  + std::to_string(t_all.size()));

    const std::vector<std::string> energies = {"10.59","10.60","10.2"};

    int ecount = 0;
    for (const auto& E : energies) {
        ++ecount;
        const std::string bc_path   = (fs::path(bincenter_json_dir)   / ("bin_centered_xsec_" + E + ".json")).string();
        const std::string pred_path = (fs::path(predictions_json_dir) / ("model_predictions_" + E + ".json")).string();

        if (!fs::exists(pred_path)) { log_warn("Missing predictions for E=" + E + " at " + pred_path + " (skipping)."); continue; }
        if (!fs::exists(bc_path))   { log_warn("Missing bin-centered data for E=" + E + " at " + bc_path + " (skipping)."); continue; }

        log_info("Energy " + E + " (" + std::to_string(ecount) + "/" + std::to_string(energies.size()) + "):");
        log_info("  loading data  " + bc_path);
        log_info("  loading preds " + pred_path);

        json j_bc   = load_json(bc_path);
        json j_pred = load_json(pred_path);
        if (j_bc.empty()   || !j_bc.contains("bins"))  { log_warn("Malformed bin-centered JSON for E=" + E + ". Skipping."); continue; }
        if (j_pred.empty() || !j_pred.contains("bins")){ log_warn("Malformed predictions JSON for E=" + E + ". Skipping."); continue; }

        for (int ix = 0; ix < (int)xB_bins.size(); ++ix) {
            const auto xb = xB_bins[ix];

            std::set<std::pair<double,double>> qs, ts;
            for (const auto& b : binning_scheme) {
                if (std::make_pair(b.xBmin,b.xBmax) == xb) {
                    qs.emplace(b.Q2min,b.Q2max);
                    ts.emplace(b.tmin,b.tmax);
                }
            }
            std::vector<std::pair<double,double>> Q2_slice(qs.begin(), qs.end());
            std::vector<std::pair<double,double>> t_slice(ts.begin(),  ts.end());
            if (Q2_slice.empty() || t_slice.empty()) continue;

            const int nrows = (int)t_slice.size();
            const int ncols = (int)Q2_slice.size();

            const int W = 320*ncols + 220;     // give more horizontal room
            const int H = 260*nrows + 260;     // extra space for title+legend

            std::ostringstream cname; cname << "c_models_"<<E<<"_xB"<<ix;
            log_info("xB slice " + std::to_string(ix) + ": grid " + std::to_string(nrows) + "x" + std::to_string(ncols));

            auto t_canvas_start = std::chrono::steady_clock::now();

            TCanvas* c = new TCanvas(cname.str().c_str(), cname.str().c_str(), W, H);
            c->cd();

            // Top band: title + single legend
            TPad* pTop = new TPad("pTop","pTop", 0.0, 0.90, 1.0, 1.0);
            pTop->SetFillStyle(0); pTop->SetBorderSize(0); pTop->Draw();
            pTop->cd();

            // Title
            TLatex head; head.SetNDC(); head.SetTextAlign(22); head.SetTextFont(42); head.SetTextSize(0.40);
            std::ostringstream tit;
            tit << "models vs. bin-centered d#sigma/d#phi"
                << Form("   %s GeV,  x_{B} #in [%.2g, %.2g]", E.c_str(), xb.first, xb.second);
            head.DrawLatex(0.50, 0.68, tit.str().c_str());

            // Global legend (compact, single row)
            // Build tiny dummy objects so the legend shows both markers and lines
            TGraph dp; dp.SetMarkerStyle(20); dp.SetMarkerColor(kBlue+1);
            TGraph dm; dm.SetMarkerStyle(25); dm.SetMarkerColor(kRed+1);
            TGraph km; km.SetLineColor(kBlue+1); km.SetLineStyle(2); km.SetLineWidth(3);
            TGraph kn; kn.SetLineColor(kRed+1);  kn.SetLineStyle(2); kn.SetLineWidth(3);
            TGraph vg; vg.SetLineColor(kBlue+1); vg.SetLineStyle(3); vg.SetLineWidth(3);
            TGraph vn; vn.SetLineColor(kRed+1);  vn.SetLineStyle(3); vn.SetLineWidth(3);
            TGraph bh; bh.SetLineColor(kBlack);  bh.SetLineStyle(2); bh.SetLineWidth(3);

            TLegend* legTop = new TLegend(0.08, 0.10, 0.92, 0.56);
            legTop->SetNColumns(5);
            legTop->SetBorderSize(0);
            legTop->SetFillStyle(0);
            legTop->SetTextFont(42);
            legTop->SetTextSize(0.22);
            legTop->AddEntry(&dp, "+ helicity (data)", "p");
            legTop->AddEntry(&dm, "- helicity (data)", "p");
            legTop->AddEntry(&km, "KM15 + (dashed)",  "l");
            legTop->AddEntry(&kn, "KM15 - (dashed)",  "l");
            legTop->AddEntry(&vg, "VGG + (dotted)",   "l");
            legTop->AddEntry(&vn, "VGG - (dotted)",   "l");
            legTop->AddEntry(&bh, "BH exact",         "l");
            legTop->Draw();

            // Grid pad below
            c->cd();
            TPad* pGrid = new TPad("pGrid","pGrid", 0.0, 0.00, 1.0, 0.90);
            pGrid->SetFillStyle(0); pGrid->SetBorderSize(0); pGrid->Draw();
            pGrid->cd();
            pGrid->Divide(ncols, nrows, 0.00, 0.00);

            // Draw panels
            for (int r=0; r<nrows; ++r) {
                const int it_global = findIndex(t_slice[r], t_all);
                if (it_global < 0) continue;
                for (int cc=0; cc<ncols; ++cc) {
                    const int iQ_global = findIndex(Q2_slice[cc], Q2_all);
                    if (iQ_global < 0) continue;

                    const std::string bkey =
                        "(" + std::to_string(ix) + "," + std::to_string(iQ_global) + "," + std::to_string(it_global) + ")";
                    log_info("Plotting panel " + bkey + "...");

                    pGrid->cd(r*ncols + cc + 1);
                    gPad->SetGrid(1,1);
                    gPad->SetTicks(1,1);
                    gPad->SetTopMargin(0.06);
                    gPad->SetBottomMargin(0.18);
                    gPad->SetLeftMargin(0.18);     // more left padding so labels do not clip
                    gPad->SetRightMargin(0.08);
                    gPad->SetLogy();

                    const json* jb_cell   = (j_bc["bins"].contains(bkey)   ? &j_bc["bins"][bkey]   : nullptr);
                    const json* jp_cell   = (j_pred["bins"].contains(bkey) ? &j_pred["bins"][bkey] : nullptr);

                    // Dynamic y-range per panel (log-safe)
                    auto [ymin_panel, ymax_panel] = panelYRangeLogSafe(jb_cell ? *jb_cell : json::object(), jp_cell);

                    TH1* frame = gPad->DrawFrame(0.0, ymin_panel, 360.0, ymax_panel);
                    frame->GetXaxis()->SetTitle("#phi (deg)");
                    frame->GetYaxis()->SetTitle("d#sigma/d#phi (bin-centered)");
                    frame->GetXaxis()->CenterTitle();
                    frame->GetYaxis()->CenterTitle();
                    frame->GetXaxis()->SetNdivisions(505);
                    frame->GetXaxis()->SetTitleSize(0.060);
                    frame->GetYaxis()->SetTitleSize(0.060);
                    frame->GetYaxis()->SetLabelSize(0.050);
                    frame->GetXaxis()->SetLabelSize(0.050);
                    frame->GetXaxis()->SetTitleOffset(1.10);
                    frame->GetYaxis()->SetTitleOffset(1.55);

                    drawDegreeTicks(gPad->GetUxmin(), gPad->GetUymin(), gPad->GetUxmax(), 0.050);

                    // Data
                    double ymax_seen = 0.0;
                    if (jb_cell){
                        if (jb_cell->contains("helicity_plus"))
                            draw_data_hel((*jb_cell)["helicity_plus"],  20, kBlue+1, frame, ymax_seen);
                        if (jb_cell->contains("helicity_minus"))
                            draw_data_hel((*jb_cell)["helicity_minus"], 25, kRed+1,  frame, ymax_seen);
                    } else {
                        log_warn("No data for " + bkey + " (panel will show models only).");
                    }

                    // Model curves
                    auto makeLine = [&](const std::vector<double>& xv, const std::vector<double>& yv,
                                        int color, int lstyle)->TGraph* {
                        auto* gr = new TGraph((int)xv.size(),
                                              const_cast<double*>(xv.data()),
                                              const_cast<double*>(yv.data()));
                        gr->SetLineColor(color);
                        gr->SetLineWidth(3);
                        gr->SetLineStyle(lstyle);
                        gr->Draw("L SAME");
                        return gr;
                    };

                    if (jp_cell){
                        std::vector<double> ph, v;
                        ph.reserve((*jp_cell)["phi_deg"].size());
                        for (size_t i=0;i<(*jp_cell)["phi_deg"].size();++i)
                            ph.push_back(jgetd((*jp_cell)["phi_deg"], i, 0.0));

                        auto pull = [&](const char* key, int color, int style)->TGraph* {
                            if (!jp_cell->contains(key)) return nullptr;
                            v.clear(); v.reserve(ph.size());
                            for (size_t i=0;i<ph.size();++i) v.push_back(jgetd((*jp_cell)[key], i, 0.0));
                            return makeLine(ph, v, color, style);
                        };

                        pull("km15_plus",  kBlue+1, 2);
                        pull("km15_minus", kRed+1,  2);
                        pull("vgg_plus",   kBlue+1, 3);
                        pull("vgg_minus",  kRed+1,  3);
                        pull("bh",         kBlack,  2);
                    } else {
                        log_warn("No predictions for " + bkey + " (panel will show data only).");
                    }

                    // per-panel label
                    TLatex lab; lab.SetNDC(); lab.SetTextSize(0.070); lab.SetTextAlign(11); lab.SetTextFont(42);
                    lab.DrawLatex(0.15, 0.94,
                        Form("Q^{2} #in [%.2g, %.2g],  -t #in [%.2g, %.2g]",
                             Q2_slice[cc].first, Q2_slice[cc].second,
                             t_slice[r].first,   t_slice[r].second));

                    log_info("Panel " + bkey + " done.");
                }
            }

            const std::string outP =
                (fs::path(output_dir)/"plots"/(std::string("models_vs_bincentered_")+E+"_xB_"+std::to_string(ix)+".png")).string();
            c->SaveAs(outP.c_str());
            delete c;

            auto t_canvas_end = std::chrono::steady_clock::now();
            auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(t_canvas_end - t_canvas_start).count();
            log_info("Saved " + outP + " (" + std::to_string(ms) + " ms).");
        }
    }

    log_info("Finished models vs. bin-centered plotting.");
}