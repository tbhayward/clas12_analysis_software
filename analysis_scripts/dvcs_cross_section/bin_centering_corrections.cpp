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
// helper: compute fbin for one (xB,Q2,t,phi-bin,helicity)
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
    // outputs
    double& fbin_km15, double& fbin_vgg
) {
    fbin_km15 = 0.0;
    fbin_vgg  = 0.0;

    // centers
    const double xB_c   = midpoint(xBmin, xBmax);
    const double Q2_c   = midpoint(Q2min, Q2max);
    const double tpos_c = midpoint(tmin_pos, tmax_pos);

    // model values at center
    const double km15_center = km15_xs(xB_c, Q2_c, tpos_c, phi_center_deg, Ebeam, hel, paths);
    const double vgg_center  =  vgg_xs(xB_c, Q2_c, tpos_c, phi_center_deg, Ebeam, hel, paths, vgg_globalfit);

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

    // accumulate over sub-bins
    double sum_km15 = 0.0; int cnt_km15 = 0;
    double sum_vgg  = 0.0; int cnt_vgg  = 0;

    for (double xb : xs){
        for (double q2 : Qs){
            for (double tp : ts){
                for (double ph : ps){
                    const double vk = km15_xs(xb, q2, tp, ph, Ebeam, hel, paths);
                    if (finite_pos(vk)) { sum_km15 += vk; ++cnt_km15; } //endif
                    const double vv =  vgg_xs(xb, q2, tp, ph, Ebeam, hel, paths, vgg_globalfit);
                    if (finite_pos(vv)) { sum_vgg  += vv; ++cnt_vgg; } //endif
                } //endfor
            } //endfor
        } //endfor
    } //endfor

    // averages
    const double avg_km15 = (cnt_km15>0 ? (sum_km15 / double(cnt_km15)) : 0.0);
    const double avg_vgg  = (cnt_vgg >0 ? (sum_vgg  / double(cnt_vgg )) : 0.0);

    if (finite_pos(avg_km15) && finite_pos(km15_center)) fbin_km15 = km15_center / avg_km15; //endif
    if (finite_pos(avg_vgg ) && finite_pos(vgg_center )) fbin_vgg  = vgg_center  / avg_vgg;  //endif
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
    bool vgg_globalfit
) {
    fs::create_directories(output_dir);
    fs::create_directories(fs::path(output_dir)/"jsons");
    fs::create_directories(fs::path(output_dir)/"plots");

    const auto xB_bins = uniqueRanges(binning_scheme, 'x');
    const auto Q2_all  = uniqueRanges(binning_scheme, 'Q');
    const auto t_all   = uniqueRanges(binning_scheme, 't');
    const auto PHI_DEG = phiCentersDeg();

    // match energies used by the RC stage
    const std::vector<std::string> energies = {"10.59","10.60","10.2"};

    for (const auto& E : energies) {
        const std::string rc_path = (fs::path(radcorr_xsec_json_dir) / ("rad_corrected_xsec_" + E + ".json")).string();
        json j_rc = load_json(rc_path);
        if (j_rc.empty() || !j_rc.contains("bins")) {
            std::cerr << "[bincenter] NOTE: missing or malformed RC file for E="<<E<<" -> "<<rc_path<<"\n";
            continue;
        } //endif

        // output container
        json jout;
        jout["energy"] = E;
        jout["bins"]   = json::object();

        // Loop bins
        for (int ix = 0; ix < (int)xB_bins.size(); ++ix) {
            const auto xb = xB_bins[ix];

            // collect Q2,t slices present in the scheme for this xB
            std::set<std::pair<double,double>> qs, ts;
            for (const auto& b : binning_scheme) {
                if (std::make_pair(b.xBmin,b.xBmax) == xb) {
                    qs.emplace(b.Q2min,b.Q2max);
                    ts.emplace(b.tmin,b.tmax);
                } //endif
            } //endfor
            std::vector<std::pair<double,double>> Q2_slice(qs.begin(), qs.end());
            std::vector<std::pair<double,double>> t_slice(ts.begin(),  ts.end());

            for (const auto& q2r : Q2_slice) {
                const int iQ_global = findIndex(q2r, Q2_all);
                if (iQ_global < 0) continue;
                for (const auto& tr : t_slice) {
                    const int it_global = findIndex(tr, t_all);
                    if (it_global < 0) continue;

                    const std::string bkey = "(" + std::to_string(ix) + "," + std::to_string(iQ_global) + "," + std::to_string(it_global) + ")";
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
                                fkm, fvg
                            );

                            // choose final fbin
                            double f_final = 1.0;
                            double f_sys   = 0.0;
                            if (finite_pos(fkm)) {
                                f_final = 0.5*(fkm + fvg);
                                f_sys   = std::fabs(fkm - fvg);
                            } else if (finite_pos(fvg)) {
                                f_final = fvg;
                                f_sys   = 1.0;
                            } else {
                                f_final = 1.0;
                                f_sys   = 0.0;
                            } //endif

                            const double x  = jgetd(xs, ip, 0.0);
                            const double ex = std::max(0.0, jgetd(xe, ip, 0.0));

                            const double yc = f_final * x;
                            const double ec = f_final * ex;  // stat error scales with factor

                            phi_bc[ip]   = phi_c;
                            y[ip]        = yc;
                            e[ip]        = ec;
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

                    jout["bins"][bkey] = {
                        {"helicity_plus",  hp_bc},
                        {"helicity_minus", hm_bc}
                    };
                } //endfor
            } //endfor
        } //endfor

        // save JSON
        const std::string out_json = (fs::path(output_dir)/"jsons"/("bin_centered_xsec_"+E+".json")).string();
        save_json(jout, out_json);
        std::cout << "[bincenter] Wrote " << out_json << "\n";

        // ---------------- PLOTS ----------------
        auto makeCanvas = [&](bool overlay){
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

                        const std::string bkey = "(" + std::to_string(ix) + "," + std::to_string(iQ_global) + "," + std::to_string(it_global) + ")";
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

                        TLatex lab;
                        lab.SetNDC(); lab.SetTextSize(0.040); lab.SetTextAlign(11);
                        lab.SetTextFont(42);
                        lab.DrawLatex(0.15, 0.94,
                            Form("Q^{2} in [%.2g, %.2g],   -t in [%.2g, %.2g]",
                                 Q2_slice[cc].first, Q2_slice[cc].second,
                                 t_slice[r].first,   t_slice[r].second));

                        if (gcp || gcm || gup || gum){
                            double x1=0.48, y1=0.70, x2=0.90, y2=0.93;
                            TLegend* leg = new TLegend(x1,y1,x2,y2);
                            leg->SetBorderSize(1);
                            leg->SetLineColor(kBlack);
                            leg->SetFillColor(kWhite);
                            leg->SetFillStyle(1001);
                            leg->SetTextFont(42);
                            leg->SetTextSize(0.038);
                            if (gcp) leg->AddEntry(gcp, "bin-centered  + helicity", "lep");
                            if (gcm) leg->AddEntry(gcm, "bin-centered  - helicity", "lep");
                            if (gup) leg->AddEntry(gup, "before BC     + helicity", "lep");
                            if (gum) leg->AddEntry(gum, "before BC     - helicity", "lep");
                            leg->Draw();
                        } //endif
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
    } //endfor

    std::cout << "[bincenter] Finished bin-centering correction.\n";
} //enddef