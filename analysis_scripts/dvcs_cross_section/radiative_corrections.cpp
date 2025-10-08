// radiative_corrections.cpp
#include "radiative_corrections.h"
#include "load_binning_scheme.h"

#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TTree.h>
#include <TGaxis.h>
#include <TH1.h>
#include <TGraphErrors.h>
#include <TPad.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace {

// ========== small style bootstrap ==========
struct StyleInit {
    StyleInit() {
        gStyle->SetOptTitle(0);
        gStyle->SetOptStat(0);
        gStyle->SetFrameLineWidth(2);
        gStyle->SetLineWidth(2);
        gStyle->SetPadTickX(1);
        gStyle->SetPadTickY(1);
        gStyle->SetLegendBorderSize(1);
        const int rf = 42;
        gStyle->SetTitleFont(rf, "XYZ");
        gStyle->SetLabelFont(rf, "XYZ");
        gStyle->SetTextFont(rf);
    }
} _style_guard;

constexpr int    N_PHI_BINS = 12;
constexpr double TWO_PI     = 2.0 * M_PI;

// --------- helpers for bin structures ----------
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
static inline std::string topoKey(const std::string& topo){
    if (topo=="(FD,FD)") return "FD_FD";
    if (topo=="(CD,FD)") return "CD_FD";
    if (topo=="(CD,FT)") return "CD_FT";
    return "FD_FD";
}
static inline std::string prettyPeriodKey(const std::string& runTagLower) {
    // runTagLower like "sp18_inb" -> "DVCS_Sp18_inb"
    std::string p = "DVCS_";
    bool upNext = true;
    for (char c : runTagLower) {
        if (upNext) { p.push_back(std::toupper(c)); upNext=false; }
        else if (c=='_') { p.push_back('_'); upNext=true; }
        else p.push_back(c);
    }
    return p;
}
static inline std::string toLower(std::string s){ std::transform(s.begin(),s.end(),s.begin(),::tolower); return s; }
static inline std::string periodToRunTagKey(const std::string& period) {
    // "DVCS_Fa18_inb" -> "fa18_inb"
    auto pos = period.find('_');
    if (pos == std::string::npos || pos+1>=period.size()) return toLower(period);
    auto tail = period.substr(pos+1);
    return toLower(tail);
}

static inline std::vector<double> phiCentersRad() {
    std::vector<double> v(N_PHI_BINS);
    const double step = TWO_PI / double(N_PHI_BINS);
    for (int i=0;i<N_PHI_BINS;++i) v[i] = (i+0.5)*step;
    return v;
}
static inline std::vector<double> phiCentersDeg() {
    std::vector<double> d(N_PHI_BINS);
    const double step = 360.0/double(N_PHI_BINS);
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

// --------- branch binder for GENERATED trees ----------
struct BranchGen {
    double x=0.0;    bool has_x=false;
    double Q2=0.0;   bool has_Q2=false;
    double t1=0.0;   bool has_t1=false; // sign as in trees
    double phi2=0.0; bool has_phi=false;
    void bind(TTree* t){
        auto bindD=[&](const char*n,double*a,bool&f){ if(t->GetBranch(n)){ t->SetBranchAddress(n,a); f=true; } };
        bindD("x",&x,has_x);
        bindD("Q2",&Q2,has_Q2);
        bindD("t1",&t1,has_t1);
        bindD("phi2",&phi2,has_phi);
    }
};

static inline int phi_index(double phi_wrapped){
    const double width = TWO_PI/double(N_PHI_BINS);
    double w = std::fmod(phi_wrapped, TWO_PI);
    if (w < 0) w += TWO_PI;
    if (w >= TWO_PI) w = std::nextafter(TWO_PI, 0.0);
    int ip = std::min(std::max(int(std::floor(w/width)),0), N_PHI_BINS-1);
    return ip;
}

// Count generated events TOT (over whole tree), and per (cell, φ)
struct GenCounts {
    double Ntot = 0.0;
    // per (ix,iQ,it) -> array[phi]
    std::map<std::tuple<int,int,int>, std::vector<double>> cell_phi;
};

static void count_generated(
    TTree* tGen,
    const std::vector<std::pair<double,double>>& xB_bins,
    const std::vector<std::pair<double,double>>& Q2_bins,
    const std::vector<std::pair<double,double>>& t_bins,
    GenCounts& out)
{
    if (!tGen) return;
    BranchGen g; g.bind(tGen);

    // prepare maps with zeros
    for (int ix=0; ix<(int)xB_bins.size(); ++ix)
    for (int iQ=0; iQ<(int)Q2_bins.size(); ++iQ)
    for (int it=0; it<(int)t_bins.size(); ++it) {
        out.cell_phi[std::make_tuple(ix,iQ,it)] = std::vector<double>(N_PHI_BINS, 0.0);
    }

    const auto findBin=[&](double v,const std::vector<std::pair<double,double>>& ranges)->int{
        for (int k=0;k<(int)ranges.size();++k) if (v>=ranges[k].first && v<ranges[k].second) return k;
        return -1;
    };

    const Long64_t n = tGen->GetEntries();
    for (Long64_t i=0;i<n;++i){
        tGen->GetEntry(i);
        if (!(g.has_x && g.has_Q2 && g.has_t1 && g.has_phi)) continue;

        out.Ntot += 1.0; // total generated events for this period/model

        const double xB = g.x;
        const double Q2 = g.Q2;
        const double tt = std::fabs(g.t1);

        const int ix = findBin(xB, xB_bins);
        const int iQ = findBin(Q2, Q2_bins);
        const int it = findBin(tt,  t_bins);

        if (ix<0||iQ<0||it<0) continue;

        const int ip = phi_index(g.phi2);
        out.cell_phi[std::make_tuple(ix,iQ,it)][ip] += 1.0;
    }
}

// --------- plotting: points with error bars, y in [0,2] ----------
static void plot_rcphi_for_period(
    const std::string& period_pretty,  // e.g. "DVCS_Sp18_inb"
    const std::vector<Binning>& binning_scheme,
    const std::vector<std::pair<double,double>>& xB_bins,
    const std::vector<std::pair<double,double>>& Q2_bins,
    const std::vector<std::pair<double,double>>& t_bins,
    const std::map<std::tuple<int,int,int>, std::vector<double>>& RC_phi,
    const std::map<std::tuple<int,int,int>, std::vector<double>>& RC_err,
    const std::string& out_dir_plots)
{
    using std::filesystem::create_directories;
    std::error_code ec;
    create_directories(out_dir_plots, ec);

    static const auto PHI_DEG = phiCentersDeg();

    for (int ix = 0; ix < (int)xB_bins.size(); ++ix) {
        const auto xb = xB_bins[ix];

        // Only Q² and |t| bins present in this x_B slice
        std::vector<std::pair<double,double>> Q2_slice, t_slice;
        {
            std::set<std::pair<double,double>> qs, ts;
            for (const auto& b : binning_scheme) {
                if (std::make_pair(b.xBmin,b.xBmax)==xb) {
                    qs.emplace(b.Q2min,b.Q2max);
                    ts.emplace(b.tmin,b.tmax);
                }
            }
            Q2_slice.assign(qs.begin(),qs.end());
            t_slice.assign(ts.begin(),ts.end());
        }
        if (Q2_slice.empty() || t_slice.empty()) continue;

        const int nrows = (int)t_slice.size();
        const int ncols = (int)Q2_slice.size();

        const int W = 280*ncols + 160;
        const int H = 240*nrows + 170;

        std::ostringstream cname; cname<<"c_rcphi_"<<period_pretty<<"_xB"<<ix;
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
        tit << Form("Radiative correction (generated)  %s   x_{B} #in [%.2g, %.2g]",
                    period_pretty.c_str(), xb.first, xb.second);
        head.DrawLatex(0.5, 0.55, tit.str().c_str());

        for (int r = 0; r < nrows; ++r) {
            const int it_global = findIndex(t_slice[r], t_bins);
            if (it_global < 0) continue;

            for (int ccol = 0; ccol < ncols; ++ccol) {
                const int iQ_global = findIndex(Q2_slice[ccol], Q2_bins);
                if (iQ_global < 0) continue;

                pGrid->cd(r*ncols + ccol + 1);
                gPad->SetGrid(1,1);
                gPad->SetTopMargin(0.08);
                gPad->SetBottomMargin(0.18);
                gPad->SetLeftMargin(0.15);
                gPad->SetRightMargin(0.10);

                // Frame: x 0..360 (consistent), y 0..2
                TH1* frame = gPad->DrawFrame(0.0, 0.0, 360.0, 2.0);
                TAxis* ax = frame->GetXaxis();
                TAxis* ay = frame->GetYaxis();

                ax->SetLabelSize(0.0001); // hide default tick labels
                ax->SetTitle("#phi (deg)");
                ay->SetTitle("rad/Born");
                ax->CenterTitle(); ay->CenterTitle();
                ax->SetNdivisions(505);
                ax->SetTitleSize(0.060);
                ay->SetTitleSize(0.060);
                ay->SetLabelSize(0.048);
                ax->SetTitleOffset(1.25);
                ay->SetTitleOffset(1.35);
                drawDegreeTicks(gPad->GetUxmin(), gPad->GetUymin(), gPad->GetUxmax(), 0.048);

                auto key = std::make_tuple(ix, iQ_global, it_global);
                auto itR  = RC_phi.find(key);
                auto itRe = RC_err.find(key);
                if (itR == RC_phi.end() || itRe == RC_err.end()) continue;

                const auto& rc  = itR->second;
                const auto& er  = itRe->second;
                const auto  PHI = phiCentersDeg();

                std::vector<double> xs, ys, eys;
                xs.reserve(rc.size()); ys.reserve(rc.size()); eys.reserve(rc.size());
                for (int ip=0; ip<N_PHI_BINS; ++ip){
                    double val = rc[ip];
                    double err = er[ip];
                    if (!std::isfinite(val)) continue;
                    xs.push_back(PHI[ip]);
                    ys.push_back(std::clamp(val, 0.0, 2.0));
                    eys.push_back((std::isfinite(err)? err : 0.0));
                }
                if (!xs.empty()){
                    TGraphErrors* gr = new TGraphErrors((int)xs.size(), xs.data(), ys.data(), nullptr, eys.data());
                    gr->SetMarkerStyle(20);
                    gr->SetMarkerSize(1.0);
                    gr->SetLineWidth(2);
                    gr->Draw("P SAME");
                }

                // panel label
                TLatex lab;
                lab.SetNDC(); lab.SetTextSize(0.040); lab.SetTextAlign(11);
                lab.SetTextFont(42);
                lab.DrawLatex(0.15, 0.94,
                    Form("Q^{2} #in [%.2g, %.2g],   -t #in [%.2g, %.2g]",
                         Q2_slice[ccol].first, Q2_slice[ccol].second,
                         t_slice[r].first,    t_slice[r].second));
            }
        }

        std::ostringstream fout;
        fout << out_dir_plots << "/rcphi_" << period_pretty << "_xB_" << ix << ".png";
        c->SaveAs(fout.str().c_str());
        delete c;
    }
}

// --------- JSON writer (per-φ arrays + raw counts & totals) ----------
static void write_period_json(
    const std::string& out_path,
    const std::vector<std::pair<double,double>>& xB_bins,
    const std::vector<std::pair<double,double>>& Q2_bins,
    const std::vector<std::pair<double,double>>& t_bins,
    double Ntot_born, double Ntot_rad,
    const std::map<std::tuple<int,int,int>, std::vector<double>>& phi_centers,
    const std::map<std::tuple<int,int,int>, std::vector<double>>& rc_vals,
    const std::map<std::tuple<int,int,int>, std::vector<double>>& rc_errs,
    const std::map<std::tuple<int,int,int>, std::vector<double>>& born_Ngen,
    const std::map<std::tuple<int,int,int>, std::vector<double>>& rad_Ngen
){
    std::ofstream ofs(out_path);
    if (!ofs) { std::cerr<<"[radcorr] Cannot open "<<out_path<<"\n"; return; }
    ofs<<std::fixed<<std::setprecision(8);
    ofs<<"{\n";
    ofs<<"  \"phi_bins\": "<<N_PHI_BINS<<",\n";
    ofs<<"  \"Ntot_born\": "<<Ntot_born<<",\n";
    ofs<<"  \"Ntot_rad\": "<<Ntot_rad<<",\n";
    ofs<<"  \"binning_meta\": {\"xB_bins\": "<<xB_bins.size()<<", \"Q2_bins\": "<<Q2_bins.size()<<", \"t_bins\": "<<t_bins.size()<<"},\n";
    ofs<<"  \"bins\": {\n";
    bool first=true;
    for (const auto& kv : rc_vals){
        if (!first) ofs<<",\n"; first=false;
        int ix,iQ,it; std::tie(ix,iQ,it)=kv.first;
        ofs<<"    \"("<<ix<<","<<iQ<<","<<it<<")\": {";
        auto arrDump=[&](const std::vector<double>& v){
            ofs<<"["; for (size_t i=0;i<v.size();++i){ if(i) ofs<<","; ofs<<v[i]; } ofs<<"]";
        };
        ofs<<"\"phi\": "; arrDump(phi_centers.at(kv.first)); ofs<<", ";
        ofs<<"\"RC\": ";  arrDump(kv.second); ofs<<", ";
        ofs<<"\"RC_err\": "; arrDump(rc_errs.at(kv.first)); ofs<<", ";
        ofs<<"\"born\": {\"Ngen\": "; arrDump(born_Ngen.at(kv.first)); ofs<<"}, ";
        ofs<<"\"rad\":  {\"Ngen\": "; arrDump(rad_Ngen.at(kv.first));  ofs<<"}";
        ofs<<"}";
    }
    ofs<<"\n  }\n}\n";
}

static void write_all_periods_json(
    const std::string& out_path,
    const std::map<std::string, std::map<std::tuple<int,int,int>, std::vector<double>>>& perPeriodRC,
    const std::map<std::string, std::map<std::tuple<int,int,int>, std::vector<double>>>& perPeriodRCErr)
{
    std::ofstream ofs(out_path);
    if (!ofs) { std::cerr<<"[radcorr] Cannot open "<<out_path<<"\n"; return; }
    ofs<<std::fixed<<std::setprecision(8);
    ofs<<"{\n  \"phi_bins\": "<<N_PHI_BINS<<",\n  \"periods\": {\n";
    bool firstP=true;
    for (const auto& pk : perPeriodRC){
        if (!firstP) ofs<<",\n"; firstP=false;
        ofs<<"    \""<<pk.first<<"\": {\"bins\":{";
        bool fb=true;
        for (const auto& kv : pk.second){
            if (!fb) ofs<<",";
            fb=false;
            int ix,iQ,it; std::tie(ix,iQ,it)=kv.first;
            ofs<<"\"("<<ix<<","<<iQ<<","<<it<<")\": {\"RC\": [";
            const auto& v = kv.second;
            for (size_t i=0;i<v.size();++i){ if(i) ofs<<","; ofs<<v[i]; }
            ofs<<"], \"RC_err\": [";
            const auto& e = perPeriodRCErr.at(pk.first).at(kv.first);
            for (size_t i=0;i<e.size(); ++i){ if(i) ofs<<","; ofs<<e[i]; }
            ofs<<"]}";
        }
        ofs<<"}}";
    }
    ofs<<"\n  }\n}\n";
}

} // anon

// =====================================================================
// Public driver
// =====================================================================
void compute_radiative_corrections(
    const std::vector<std::string>& periods,
    const std::vector<std::string>& topologies,                     // kept for API compat; UNUSED here
    const std::vector<Binning>& binning_scheme,
    const std::map<std::string, TTree*>& genMcTrees_norad,          // keys: "sp18_inb_gen", ...
    const std::map<std::string, TTree*>& recMcTrees_norad,          // UNUSED
    const std::map<std::string, TTree*>& genMcTrees_rad,            // keys: "sp18_inb_gen_rad", ...
    const std::map<std::string, TTree*>& recMcTrees_rad,            // UNUSED
    const std::string& combined_cuts_json_path,                     // UNUSED
    const std::string& out_root_dir)
{
    (void)topologies;
    (void)recMcTrees_norad;
    (void)recMcTrees_rad;
    (void)combined_cuts_json_path;

    namespace fs = std::filesystem;

    const auto xB_bins = uniqueRanges(binning_scheme, 'x');
    const auto Q2_bins = uniqueRanges(binning_scheme, 'Q');
    const auto t_bins  = uniqueRanges(binning_scheme, 't');

    // Output directories
    const fs::path json_dir  = fs::path(out_root_dir)/"jsons";
    const fs::path plot_root = fs::path(out_root_dir)/"radiative_correction_plots";
    std::error_code ec;
    fs::create_directories(json_dir, ec);
    fs::create_directories(plot_root, ec);

    // For all-periods summary
    std::map<std::string, std::map<std::tuple<int,int,int>, std::vector<double>>> perPeriodRC;
    std::map<std::string, std::map<std::tuple<int,int,int>, std::vector<double>>> perPeriodRCErr;

    const auto PHI_RAD = phiCentersRad();

    for (const auto& period : periods) {
        const std::string runTag = periodToRunTagKey(period);            // "fa18_inb", "fa18_inb_supp", ...
        const std::string pretty = prettyPeriodKey(runTag);              // "DVCS_Fa18_inb", etc.
        const std::string baseRunTag = (runTag=="fa18_inb_supp") ? "fa18_inb" : runTag;

        auto getOrNull = [](const auto& m, const std::string& k)->TTree*{
            auto it=m.find(k); return (it!=m.end()? it->second : nullptr);
        };

        TTree* gB = getOrNull(genMcTrees_norad, baseRunTag + "_gen");
        TTree* gR = getOrNull(genMcTrees_rad,   baseRunTag + "_gen_rad");

        if (!gB || !gR) {
            std::cerr<<"[radcorr] Missing generated MC trees (born or rad) for "<<baseRunTag<<"\n";
            continue;
        }

        // Count totals and per-cell,per-phi
        GenCounts cntB, cntR;
        count_generated(gB, xB_bins, Q2_bins, t_bins, cntB);
        count_generated(gR, xB_bins, Q2_bins, t_bins, cntR);

        // Per-cell arrays for JSON & plots
        std::map<std::tuple<int,int,int>, std::vector<double>> phi_centers_map;
        std::map<std::tuple<int,int,int>, std::vector<double>> RC_phi_map, RC_err_map;
        std::map<std::tuple<int,int,int>, std::vector<double>> bornNgen_map, radNgen_map;

        // Build RC(phi) = (a/NrTot) / (b/NbTot), with multinomial proportion errors
        for (int ix=0; ix<(int)xB_bins.size(); ++ix)
        for (int iQ=0; iQ<(int)Q2_bins.size(); ++iQ)
        for (int it=0; it<(int)t_bins.size(); ++it) {
            auto key3 = std::make_tuple(ix,iQ,it);

            const auto& nb = cntB.cell_phi.at(key3);
            const auto& nr = cntR.cell_phi.at(key3);

            std::vector<double> RC (N_PHI_BINS, std::numeric_limits<double>::quiet_NaN());
            std::vector<double> eRC(N_PHI_BINS, std::numeric_limits<double>::quiet_NaN());

            for (int ip=0; ip<N_PHI_BINS; ++ip){
                const double a = nr[ip]; // rad in cell, phi
                const double b = nb[ip]; // born in cell, phi
                const double NrTot = std::max(cntR.Ntot, 0.0);
                const double NbTot = std::max(cntB.Ntot, 0.0);

                if (NrTot <= 0.0 || NbTot <= 0.0 || b <= 0.0) {
                    RC[ip]  = std::numeric_limits<double>::quiet_NaN();
                    eRC[ip] = std::numeric_limits<double>::quiet_NaN();
                    continue;
                }

                const double pr_hat = a / NrTot;
                const double pb_hat = b / NbTot;

                if (pb_hat <= 0.0) {
                    RC[ip]  = std::numeric_limits<double>::quiet_NaN();
                    eRC[ip] = std::numeric_limits<double>::quiet_NaN();
                    continue;
                }

                const double R = pr_hat / pb_hat;

                // multinomial proportion variance approx: var(p̂)=p̂(1-p̂)/N
                const double var_pr = pr_hat * std::max(0.0, 1.0 - pr_hat) / std::max(NrTot, 1.0);
                const double var_pb = pb_hat * std::max(0.0, 1.0 - pb_hat) / std::max(NbTot, 1.0);

                double sR = std::numeric_limits<double>::quiet_NaN();
                if (pr_hat>0.0 && pb_hat>0.0) {
                    const double rel2 = (var_pr/(pr_hat*pr_hat)) + (var_pb/(pb_hat*pb_hat));
                    sR = std::abs(R) * std::sqrt(std::max(0.0, rel2));
                }

                RC[ip]  = R;
                eRC[ip] = sR;
            }

            phi_centers_map[key3] = PHI_RAD;
            RC_phi_map[key3]      = RC;
            RC_err_map[key3]      = eRC;
            bornNgen_map[key3]    = cntB.cell_phi.at(key3);
            radNgen_map[key3]     = cntR.cell_phi.at(key3);
        }

        // Save per-period JSON
        {
            const fs::path outP = fs::path(out_root_dir)/"jsons"/("radiative_corrections_"+pretty+".json");
            write_period_json(outP.string(), xB_bins, Q2_bins, t_bins,
                              cntB.Ntot, cntR.Ntot,
                              phi_centers_map, RC_phi_map, RC_err_map,
                              bornNgen_map,    radNgen_map);
            std::cout<<"[radcorr] Wrote "<<outP.string()<<"\n";
        }

        // Plots
        const fs::path plots_dir = fs::path(out_root_dir)/"radiative_correction_plots"/runTag;
        fs::create_directories(plots_dir, ec);
        plot_rcphi_for_period(pretty, binning_scheme, xB_bins, Q2_bins, t_bins,
                              RC_phi_map, RC_err_map, plots_dir.string());

        // All-periods summary stash
        perPeriodRC[pretty]    = std::move(RC_phi_map);
        perPeriodRCErr[pretty] = std::move(RC_err_map);
    }

    // all periods file
    write_all_periods_json((fs::path(out_root_dir)/"jsons"/"radiative_corrections_all_periods.json").string(),
                           perPeriodRC, perPeriodRCErr);

    std::cout<<"[radcorr] Radiative corrections (generated-only, per-phi) complete.\n";
}