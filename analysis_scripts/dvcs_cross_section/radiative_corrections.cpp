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
#include <TGraph.h>
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
static inline std::string periodToRunTagKey(const std::string& period) {
    // "DVCS_Fa18_inb" -> "fa18_inb"
    auto pos = period.find('_');
    if (pos == std::string::npos || pos+1>=period.size()) {
        std::string s=period; std::transform(s.begin(),s.end(),s.begin(),::tolower); return s;
    }
    std::string tail = period.substr(pos+1);
    std::transform(tail.begin(), tail.end(), tail.begin(), ::tolower);
    return tail;
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

// --------- exclusivity cut structures ----------
struct Stats { double mean=0.0, std=0.0; };
struct CutDict { std::map<std::string, Stats> data; std::map<std::string, Stats> mc; };

// Very light JSON puller to get the MC means/stds we need from combined_cuts.json
// Expect keys like: "DVCS_Sp18_inb_FD_FD": {"data":{...},"mc":{...}}
static bool load_combined_cuts_mc(const std::string& path,
                                  std::map<std::string, CutDict>& out)
{
    std::ifstream ifs(path);
    if (!ifs) { std::cerr<<"[radcorr] Cannot open cuts json "<<path<<"\n"; return false; }
    std::string s((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());

    size_t pos=0;
    while (true){
        size_t q1 = s.find('"', pos); if (q1==std::string::npos) break;
        size_t q2 = s.find('"', q1+1); if (q2==std::string::npos) break;
        std::string key = s.substr(q1+1, q2-q1-1); // e.g. DVCS_Sp18_inb_FD_FD
        size_t objS = s.find('{', q2); if (objS==std::string::npos) break;
        int d=0; size_t j=objS;
        for (; j<s.size(); ++j){ if(s[j]=='{') d++; else if(s[j]=='}'){ d--; if(!d){ ++j; break; } } }
        std::string obj = s.substr(objS, j-objS);

        auto parseStats = [&](const std::string& blk)->std::map<std::string,Stats>{
            std::map<std::string,Stats> m;
            size_t p=0;
            while (true){
                size_t k1 = blk.find('"', p); if (k1==std::string::npos) break;
                size_t k2 = blk.find('"', k1+1); if (k2==std::string::npos) break;
                std::string vname = blk.substr(k1+1, k2-k1-1);
                size_t vs = blk.find('{', k2); if (vs==std::string::npos) break;
                int d2=0; size_t jj=vs;
                for (; jj<blk.size(); ++jj){ if(blk[jj]=='{') d2++; else if(blk[jj]=='}'){ d2--; if(!d2){ ++jj; break; } } }
                std::string sobj = blk.substr(vs, jj-vs);
                auto findD=[&](const char* pat)->double{
                    size_t pp=sobj.find(pat); if (pp==std::string::npos) return 0.0;
                    pp=sobj.find(':',pp); if (pp==std::string::npos) return 0.0;
                    size_t a=pp+1; while (a<sobj.size() && isspace((unsigned char)sobj[a])) ++a;
                    size_t b=a; while (b<sobj.size() &&
                        (isdigit((unsigned char)sobj[b])||sobj[b]=='-'||sobj[b]=='.'||sobj[b]=='e'||sobj[b]=='E'||sobj[b]=='+')) ++b;
                    try { return std::stod(sobj.substr(a,b-a)); } catch(...) { return 0.0; }
                };
                Stats st; st.mean = findD("\"mean\""); st.std = std::fabs(findD("\"std\""));
                m[vname]=st;
                p=jj;
            }
            return m;
        };

        CutDict cd;
        // find "data":{...}
        size_t pd = obj.find("\"data\""); if (pd!=std::string::npos){
            size_t ds = obj.find('{', pd);
            int dd=0; size_t jd=ds;
            for (; jd<obj.size(); ++jd){ if(obj[jd]=='{') dd++; else if(obj[jd]=='}'){ dd--; if(!dd){ ++jd; break; } } }
            if (ds!=std::string::npos) cd.data = parseStats(obj.substr(ds, jd-ds));
        }
        // find "mc":{...}
        size_t pm = obj.find("\"mc\""); if (pm!=std::string::npos){
            size_t ms = obj.find('{', pm);
            int dm=0; size_t jm=ms;
            for (; jm<obj.size(); ++jm){ if(obj[jm]=='{') dm++; else if(obj[jm]=='}'){ dm--; if(!dm){ ++jm; break; } } }
            if (ms!=std::string::npos) cd.mc = parseStats(obj.substr(ms, jm-ms));
        }

        out[key]=std::move(cd);
        pos=j;
    }
    return !out.empty();
}

static inline bool within3Sigma(double val, const Stats& s) {
    return (val >= s.mean - 3.0 * s.std) && (val <= s.mean + 3.0 * s.std);
}

// --------- branch binder (reco variables, same names as earlier) ----------
struct Branch {
    int detector1=0, detector2=0; bool has_d1=false, has_d2=false;
    double t1=0.0;                bool has_t1=false;
    double open_angle_ep2=0.0;    bool has_oa=false;
    double Emiss2=0.0;            bool has_Emiss2=false;
    double Mx2=0.0;               bool has_Mx2=false;
    double Mx2_1=0.0;             bool has_Mx2_1=false;
    double Mx2_2=0.0;             bool has_Mx2_2=false;
    double pTmiss=0.0;            bool has_pT=false;
    double xF=0.0;                bool has_xF=false;
    double Delta_phi=0.0;         bool has_Dphi=false;
    double x=0.0;                 bool has_x=false;
    double Q2=0.0;                bool has_Q2=false;
    double t1abs=0.0;             // |t|
    void bind(TTree* t){
        auto bindI=[&](const char*n,int*a,bool&f){ if(t->GetBranch(n)){ t->SetBranchAddress(n,a); f=true; } };
        auto bindD=[&](const char*n,double*a,bool&f){ if(t->GetBranch(n)){ t->SetBranchAddress(n,a); f=true; } };
        bindI("detector1",&detector1,has_d1);
        bindI("detector2",&detector2,has_d2);
        bindD("t1",&t1,has_t1);
        bindD("open_angle_ep2",&open_angle_ep2,has_oa);
        bindD("Emiss2",&Emiss2,has_Emiss2);
        bindD("Mx2",&Mx2,has_Mx2);
        bindD("Mx2_1",&Mx2_1,has_Mx2_1);
        bindD("Mx2_2",&Mx2_2,has_Mx2_2);
        bindD("pTmiss",&pTmiss,has_pT);
        bindD("xF",&xF,has_xF);
        bindD("Delta_phi",&Delta_phi,has_Dphi);
        bindD("x",&x,has_x);
        bindD("Q2",&Q2,has_Q2);
    }
    inline bool passesTopo(const std::string& topo) const {
        if (topo=="(FD,FD)") return has_d1 && has_d2 && detector1==1 && detector2==1;
        if (topo=="(CD,FD)") return has_d1 && has_d2 && detector1==2 && detector2==1;
        if (topo=="(CD,FT)") return has_d1 && has_d2 && detector1==2 && detector2==0;
        return false;
    }
    inline bool passesGlobalKin() const {
        if (!(has_oa && has_t1 && has_pT)) return false;
        if (open_angle_ep2 <= 5.0) return false;
        if ((-t1) > 1.0) return false;
        if (pTmiss > 0.20) return false;
        return true;
    }
    std::map<std::string,double> valmap() const {
        std::map<std::string,double> m;
        if (has_Dphi) m["Delta_phi"]=Delta_phi;
        if (has_oa)   m["theta_gamma_gamma"]=0.0; // not used directly; kept for completeness
        if (has_pT)   m["pTmiss"]=pTmiss;
        if (has_xF)   m["xF"]=xF;
        if (has_Emiss2) m["Emiss2"]=Emiss2;
        if (has_Mx2)    m["Mx2"]=Mx2;
        if (has_Mx2_1)  m["Mx2_1"]=Mx2_1;
        if (has_Mx2_2)  m["Mx2_2"]=Mx2_2;
        return m;
    }
};

// apply 3σ using the **MC** side of the cuts
static bool passes3SigmaMc(const std::map<std::string, Stats>& mcCuts,
                           const std::map<std::string,double>& values)
{
    for (const auto& kv : mcCuts) {
        auto it=values.find(kv.first);
        if (it==values.end()) continue;
        if (!within3Sigma(it->second, kv.second)) return false;
    }
    return true;
}

// --------- bin counting core (integrated over φ) ----------
struct Counts3D { double norad=0.0; double rad=0.0; };

static void count_in_bins_for_tree(
    TTree* t,
    const std::string& topo,
    const std::vector<std::pair<double,double>>& xB_bins,
    const std::vector<std::pair<double,double>>& Q2_bins,
    const std::vector<std::pair<double,double>>& t_bins,
    const std::map<std::string,Stats>& mcCuts,
    std::map<std::tuple<int,int,int>, double>& accum)
{
    if (!t) return;
    Branch b; b.bind(t);
    auto findBin=[&](double v,const std::vector<std::pair<double,double>>& ranges)->int{
        for (int k=0;k<(int)ranges.size();++k) if (v>=ranges[k].first && v<ranges[k].second) return k;
        return -1;
    };

    const Long64_t n = t->GetEntries();
    for (Long64_t i=0;i<n;++i){
        t->GetEntry(i);
        if (!b.passesTopo(topo)) continue;
        if (!b.passesGlobalKin()) continue;
        if (!(b.has_x && b.has_Q2 && b.has_t1)) continue;

        std::map<std::string,double> vals = b.valmap();
        if (!passes3SigmaMc(mcCuts, vals)) continue;

        const double xB = b.x;
        const double Q2 = b.Q2;
        const double tt = std::fabs(b.t1);

        int ix = findBin(xB, xB_bins);
        int iQ = findBin(Q2, Q2_bins);
        int it = findBin(tt, t_bins);
        if (ix<0||iQ<0||it<0) continue;

        accum[std::make_tuple(ix,iQ,it)] += 1.0;
    }
}

// --------- plotting (one constant line per panel, y in [0,2]) ----------
static void plot_rc_for_period(
    const std::string& period_pretty,  // e.g. "DVCS_Sp18_inb"
    const std::vector<Binning>& binning_scheme,
    const std::vector<std::pair<double,double>>& xB_bins,
    const std::vector<std::pair<double,double>>& Q2_bins,
    const std::vector<std::pair<double,double>>& t_bins,
    const std::map<std::tuple<int,int,int>, double>& RC,
    const std::string& out_dir_plots)
{
    using std::filesystem::create_directories;
    std::error_code ec;
    create_directories(out_dir_plots, ec);

    static const auto PHI_DEG = phiCentersDeg();

    for (int ix = 0; ix < (int)xB_bins.size(); ++ix) {
        const auto xb = xB_bins[ix];

        // Only the Q² and |t| bins that exist in this xB slice
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

        std::ostringstream cname; cname<<"c_rc_"<<period_pretty<<"_xB"<<ix;
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
        tit << Form("Radiative correction  %s   x_{B} #in [%.2g, %.2g]",
                    period_pretty.c_str(), xb.first, xb.second);
        head.DrawLatex(0.5, 0.55, tit.str().c_str());

        // Panels
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

                // frame: x 0..360 (for consistent look), y 0..2
                TH1* frame = gPad->DrawFrame(0.0, 0.0, 360.0, 2.0);
                TAxis* ax = frame->GetXaxis();
                TAxis* ay = frame->GetYaxis();

                ax->SetLabelSize(0.0001); // hide default tick labels, overlay custom
                ax->SetTitle("#phi (deg)");
                ay->SetTitle("R_{C} = N_{rad}/N_{norad}");
                ax->CenterTitle(); ay->CenterTitle();
                ax->SetNdivisions(505);
                ax->SetTitleSize(0.060);
                ay->SetTitleSize(0.060);
                ay->SetLabelSize(0.048);
                ax->SetTitleOffset(1.25);
                ay->SetTitleOffset(1.35);

                drawDegreeTicks(gPad->GetUxmin(), gPad->GetUymin(), gPad->GetUxmax(), 0.048);

                auto itCell = RC.find(std::make_tuple(ix, iQ_global, it_global));
                if (itCell == RC.end()) continue;

                const double rc = itCell->second;

                // draw as a thin horizontal line
                const int NS = 2;
                double xs[NS] = {0.0, 360.0};
                double ys[NS] = {rc, rc};
                TGraph* g = new TGraph(NS, xs, ys);
                g->SetLineColor(kRed);
                g->SetLineWidth(2);
                g->Draw("L SAME");

                // annotation
                TLatex lab;
                lab.SetNDC(); lab.SetTextSize(0.040); lab.SetTextAlign(11);
                lab.SetTextFont(42);
                lab.DrawLatex(0.15, 0.94,
                    Form("Q^{2} #in [%.2g, %.2g],   -t #in [%.2g, %.2g]",
                         Q2_slice[ccol].first, Q2_slice[ccol].second,
                         t_slice[r].first,    t_slice[r].second));

                // legend (opaque)
                TLegend* leg = new TLegend(0.60, 0.76, 0.90, 0.92);
                leg->SetBorderSize(1);
                leg->SetLineColor(kBlack);
                leg->SetFillColor(kWhite);
                leg->SetFillStyle(1001);
                leg->SetTextFont(42);
                leg->SetTextSize(0.040);
                leg->AddEntry((TObject*)nullptr, Form("R_{C} = %.3f", rc), "");
                leg->Draw();
            }
        }

        std::ostringstream fout;
        fout << out_dir_plots << "/rc_" << period_pretty << "_xB_" << ix << ".png";
        c->SaveAs(fout.str().c_str());
        delete c;
    }
}

// --------- JSON writers ----------
static void write_period_json(
    const std::string& out_path,
    const std::vector<std::pair<double,double>>& xB_bins,
    const std::vector<std::pair<double,double>>& Q2_bins,
    const std::vector<std::pair<double,double>>& t_bins,
    const std::map<std::tuple<int,int,int>, double>& rcVals,
    const std::map<std::tuple<int,int,int>, std::pair<double,double>>& rawCounts // (norad,rad)
){
    std::ofstream ofs(out_path);
    if (!ofs) { std::cerr<<"[radcorr] Cannot open "<<out_path<<"\n"; return; }
    ofs<<std::fixed<<std::setprecision(8);
    ofs<<"{\n";
    ofs<<"  \"binning_meta\": {\"xB_bins\": "<<xB_bins.size()<<", \"Q2_bins\": "<<Q2_bins.size()<<", \"t_bins\": "<<t_bins.size()<<"},\n";
    ofs<<"  \"bins\": {\n";
    bool first=true;
    for (const auto& kv : rcVals){
        if (!first) ofs<<",\n"; first=false;
        int ix,iQ,it; std::tie(ix,iQ,it)=kv.first;
        double rc = kv.second;
        auto itc = rawCounts.find(kv.first);
        double nN=0.0, nR=0.0;
        if (itc!=rawCounts.end()){ nN=itc->second.first; nR=itc->second.second; }
        ofs<<"    \"("<<ix<<","<<iQ<<","<<it<<")\": {\"RC\": "<<rc<<", \"N_norad\": "<<nN<<", \"N_rad\": "<<nR<<"}";
    }
    ofs<<"\n  }\n}\n";
}

static void write_all_periods_json(
    const std::string& out_path,
    const std::map<std::string, std::map<std::tuple<int,int,int>, double>>& perPeriodRC)
{
    std::ofstream ofs(out_path);
    if (!ofs) { std::cerr<<"[radcorr] Cannot open "<<out_path<<"\n"; return; }
    ofs<<std::fixed<<std::setprecision(8);
    ofs<<"{\n  \"periods\": {\n";
    bool firstP=true;
    for (const auto& pk : perPeriodRC){
        if (!firstP) ofs<<",\n"; firstP=false;
        ofs<<"    \""<<pk.first<<"\": {";
        const auto& rc = pk.second;
        ofs<<"\"bins\":{";
        bool fb=true;
        for (const auto& kv : rc){
            if (!fb) ofs<<",";
            fb=false;
            int ix,iQ,it; std::tie(ix,iQ,it)=kv.first;
            ofs<<"\"("<<ix<<","<<iQ<<","<<it<<")\":"<<kv.second;
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
    const std::vector<std::string>& periods,                        // e.g. DVCS_Fa18_inb, ...
    const std::vector<std::string>& topologies,                     // {"(FD,FD)","(CD,FD)","(CD,FT)"}
    const std::vector<Binning>& binning_scheme,
    const std::map<std::string, TTree*>& recMcTrees_norad,          // keys: "sp18_inb_rec", ...
    const std::map<std::string, TTree*>& recMcTrees_rad,            // keys: "sp18_inb_rec_rad", ...
    const std::string& combined_cuts_json_path,                     // "output/jsons/combined_cuts.json"
    const std::string& out_root_dir)                                // "output"
{
    namespace fs = std::filesystem;

    const auto xB_bins = uniqueRanges(binning_scheme, 'x');
    const auto Q2_bins = uniqueRanges(binning_scheme, 'Q');
    const auto t_bins  = uniqueRanges(binning_scheme, 't');

    // Load MC-side cuts
    std::map<std::string, CutDict> cuts;
    if (!load_combined_cuts_mc(combined_cuts_json_path, cuts)) {
        std::cerr<<"[radcorr] WARNING: could not read combined cuts; proceeding without 3σ exclusivity cuts.\n";
    }

    // Output directories
    const fs::path json_dir  = fs::path(out_root_dir)/"jsons";
    const fs::path plot_root = fs::path(out_root_dir)/"radiative_correction_plots";
    std::error_code ec;
    fs::create_directories(json_dir, ec);
    fs::create_directories(plot_root, ec);

    // Storage for the "all-periods" file
    std::map<std::string, std::map<std::tuple<int,int,int>, double>> perPeriodRC;

    for (const auto& period : periods) {
        const std::string runTag = periodToRunTagKey(period);          // "fa18_inb", "fa18_inb_supp", ...
        const std::string pretty = prettyPeriodKey(runTag);            // "DVCS_Fa18_inb", etc.

        // pick which actual key to use for Fa18_inb_supp (reuse Fa18_inb)
        const std::string baseRunTag = (runTag=="fa18_inb_supp") ? "fa18_inb" : runTag;

        // accumulate counts across requested topologies, then make RC
        std::map<std::tuple<int,int,int>, double> counts_norad; // per (ix,iQ,it)
        std::map<std::tuple<int,int,int>, double> counts_rad;

        for (const auto& topo : topologies) {
            // MC-side cuts key from combined_cuts.json:
            // "DVCS_Sp18_inb_FD_FD", "DVCS_Fa18_inb_FD_FD", ...
            std::string cutsKey = prettyPeriodKey(baseRunTag) + "_" + topoKey(topo);
            std::map<std::string,Stats> mcCuts;
            if (cuts.count(cutsKey)) mcCuts = cuts[cutsKey].mc;

            auto getOrNull = [](const auto& m, const std::string& k)->TTree*{
                auto it=m.find(k); return (it!=m.end()? it->second : nullptr);
            };

            TTree* tNorad = getOrNull(recMcTrees_norad, baseRunTag + "_rec");
            TTree* tRad   = getOrNull(recMcTrees_rad,   baseRunTag + "_rec_rad");

            if (!tNorad || !tRad) {
                std::cerr<<"[radcorr] Missing MC trees (norad or rad) for "<<baseRunTag<<" topo "<<topo<<"\n";
                continue;
            }

            count_in_bins_for_tree(tNorad, topo, xB_bins, Q2_bins, t_bins, mcCuts, counts_norad);
            count_in_bins_for_tree(tRad,   topo, xB_bins, Q2_bins, t_bins, mcCuts, counts_rad);
        }

        // build RC per bin (with tiny floor to avoid 0/0)
        std::map<std::tuple<int,int,int>, double> RC;
        std::map<std::tuple<int,int,int>, std::pair<double,double>> rawCounts;
        for (int ix=0; ix<(int)xB_bins.size(); ++ix)
        for (int iQ=0; iQ<(int)Q2_bins.size(); ++iQ)
        for (int it=0; it<(int)t_bins.size(); ++it) {
            auto key = std::make_tuple(ix,iQ,it);
            double Nn = counts_norad.count(key)? counts_norad[key] : 0.0;
            double Nr = counts_rad.count(key)?   counts_rad[key]   : 0.0;

            // Handle empty bins gracefully: RC = 1.0 if both zero; else (Nr / Nn) with epsilon
            double rc = 1.0;
            if (Nn<=0.0 && Nr<=0.0) rc = 1.0;
            else if (Nn<=0.0 && Nr>0.0) rc = 2.0; // cap visually later
            else rc = Nr / std::max(Nn, 1e-12);

            // Clip to plotting range a bit so legends & downstream don't blow up
            if (!std::isfinite(rc)) rc = 1.0;
            rc = std::clamp(rc, 0.0, 2.0);

            RC[key]=rc;
            rawCounts[key] = {Nn,Nr};
        }

        // Save per-period JSON
        {
            const fs::path outP = json_dir/("radiative_corrections_"+pretty+".json");
            write_period_json(outP.string(), xB_bins, Q2_bins, t_bins, RC, rawCounts);
            std::cout<<"[radcorr] Wrote "<<outP.string()<<"\n";
        }

        // Plots directory for this period
        const fs::path plots_dir = plot_root / runTag; // keep original runTag for directory
        fs::create_directories(plots_dir, ec);
        plot_rc_for_period(pretty, binning_scheme, xB_bins, Q2_bins, t_bins, RC, plots_dir.string());

        // Stash for "all periods"
        perPeriodRC[pretty]=std::move(RC);

        // If we processed fa18_inb and the requested *period* is fa18_inb_supp, we already
        // used baseRunTag=fa18_inb. For completeness, also write the duplicate JSON/plots
        // if the *period* itself is the supplemental (directory name changes).
        if (runTag=="fa18_inb_supp") {
            // nothing extra: we already keyed plots_dir with runTag so they land under fa18_inb_supp
        }
    }

    // all periods file
    write_all_periods_json((fs::path(out_root_dir)/"jsons"/"radiative_corrections_all_periods.json").string(),
                           perPeriodRC);

    std::cout<<"[radcorr] Radiative corrections complete.\n";
}