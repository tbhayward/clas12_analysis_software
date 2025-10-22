// total_counts.cpp
// ------------------------------------------------------------
// Compute helicity-separated total counts after exclusivity cuts.
// - Reads: combined_cuts.json (3σ means/stds)
// - Reads: DVCS data trees per period
// - Writes:
//     • output/jsons/total_counts.json             (all groups, nested "groups")
//     • output/jsons/total_counts_<group>.json     (flat per-group file)
// - Produces per-group φ-binned plots under output/total_counts_plots/<group>/
//
// Each group (period) JSON contains:
// {
//   "binning_meta": { ... },
//   "bins": {
//      "(ix,iQ2,it,ip)": {
//         "helicity": { "+1": N_plus, "-1": N_minus },
//         "total": N_total
//      }, ...
//   }
// }
//
// Combined master JSON structure:
// {
//   "binning_meta": {...},
//   "groups": { "<group>": { "bins": {...} }, ... }
// }
// ------------------------------------------------------------

#include "total_counts.h"
#include <TTree.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TPad.h>
#include <TH1.h>
#include <TStyle.h>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

// ------------------------------------------------------------
// Helpers and constants
// ------------------------------------------------------------
static constexpr int N_PHI_BINS = 12;
static constexpr double TWO_PI = 2.0 * M_PI;

static inline std::string toLower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c){ return std::tolower(c); });
    return s;
}

static inline double wrapToTwoPi(double phi) {
    if (!std::isfinite(phi)) return 0.0;
    double w = std::fmod(phi, TWO_PI);
    if (w < 0.0) w += TWO_PI;
    if (w >= TWO_PI) w = std::nextafter(TWO_PI, 0.0);
    return w;
}

static inline int phiToBin(double phi) {
    double w = wrapToTwoPi(phi);
    double width = TWO_PI / static_cast<double>(N_PHI_BINS);
    int idx = static_cast<int>(std::floor(w / width));
    if (idx < 0) idx = 0;
    if (idx >= N_PHI_BINS) idx = N_PHI_BINS - 1;
    return idx;
}

static std::vector<std::pair<double,double>> uniqueRanges(const std::vector<Binning>& scheme, char which) {
    std::set<std::pair<double,double>> s;
    for (const auto& b : scheme) {
        if (which == 'x') s.emplace(b.xBmin, b.xBmax);
        else if (which == 'Q') s.emplace(b.Q2min, b.Q2max);
        else if (which == 't') s.emplace(b.tmin, b.tmax);
    }
    return std::vector<std::pair<double,double>>(s.begin(), s.end());
}

static int findBin(double v, const std::vector<std::pair<double,double>>& ranges) {
    for (int i = 0; i < static_cast<int>(ranges.size()); ++i)
        if (v >= ranges[i].first && v < ranges[i].second) return i;
    return -1;
}

static inline std::string keyStr(int ix,int iQ,int it,int ip) {
    std::ostringstream os; os<<"("<<ix<<","<<iQ<<","<<it<<","<<ip<<")";
    return os.str();
}

struct HelCounts { long long plus=0, minus=0; };

// ------------------------------------------------------------
// Write one flat JSON per group (top-level "bins")
// ------------------------------------------------------------
static void write_total_counts_group_json(
    const std::string& out_path,
    const std::map<std::tuple<int,int,int,int>, HelCounts>& table,
    int nPhi,
    const std::vector<std::pair<double,double>>& xB_bins,
    const std::vector<std::pair<double,double>>& Q2_bins,
    const std::vector<std::pair<double,double>>& t_bins)
{
    std::ofstream ofs(out_path);
    if (!ofs) { std::cerr << "[total_counts][ERROR] Cannot open " << out_path << "\n"; return; }
    ofs << std::fixed << std::setprecision(8);
    ofs << "{\n";
    ofs << "  \"binning_meta\": {\"phi_bins\": " << nPhi
        << ", \"xB_bins\": " << xB_bins.size()
        << ", \"Q2_bins\": " << Q2_bins.size()
        << ", \"t_bins\": "  << t_bins.size() << "},\n";
    ofs << "  \"bins\": {\n";
    bool first=true;
    for (const auto& kv : table) {
        if (!first) ofs << ",\n"; first=false;
        const HelCounts& hc = kv.second;
        int ix,iq,itb,ip; std::tie(ix,iq,itb,ip) = kv.first;
        ofs << "    \"" << keyStr(ix,iq,itb,ip) << "\": {"
            << "\"helicity\":{\"+1\":" << hc.plus
            << ",\"-1\":" << hc.minus << "},"
            << "\"total\":" << (hc.plus + hc.minus)
            << "}";
    }
    ofs << "\n  }\n}\n";
    ofs.close();
    std::cout << "[total_counts] Wrote " << out_path << "\n";
}

// ------------------------------------------------------------
// Write master JSON (nested "groups")
// ------------------------------------------------------------
static void write_total_counts_master_json(
    const std::string& out_path,
    const std::map<std::string, std::map<std::tuple<int,int,int,int>, HelCounts>>& allGroups,
    int nPhi,
    const std::vector<std::pair<double,double>>& xB_bins,
    const std::vector<std::pair<double,double>>& Q2_bins,
    const std::vector<std::pair<double,double>>& t_bins)
{
    std::ofstream ofs(out_path);
    if (!ofs) { std::cerr << "[total_counts][ERROR] Cannot open " << out_path << "\n"; return; }
    ofs << std::fixed << std::setprecision(8);
    ofs << "{\n";
    ofs << "  \"binning_meta\": {\"phi_bins\": " << nPhi
        << ", \"xB_bins\": " << xB_bins.size()
        << ", \"Q2_bins\": " << Q2_bins.size()
        << ", \"t_bins\": "  << t_bins.size() << "},\n";
    ofs << "  \"groups\": {\n";

    bool firstG = true;
    for (const auto& gkv : allGroups) {
        if (!firstG) ofs << ",\n";
        firstG = false;
        ofs << "    \"" << gkv.first << "\": { \"bins\": {\n";
        bool firstB = true;
        for (const auto& kv : gkv.second) {
            if (!firstB) ofs << ",\n"; firstB=false;
            const HelCounts& hc = kv.second;
            int ix,iq,itb,ip; std::tie(ix,iq,itb,ip)=kv.first;
            ofs << "      \"" << keyStr(ix,iq,itb,ip) << "\": {"
                << "\"helicity\":{\"+1\":" << hc.plus
                << ",\"-1\":" << hc.minus << "},"
                << "\"total\":" << (hc.plus + hc.minus)
                << "}";
        }
        ofs << "\n    }}";
    }
    ofs << "\n  }\n}\n";
    ofs.close();
    std::cout << "[total_counts] Wrote master " << out_path << "\n";
}

// ------------------------------------------------------------
// Plotting helper
// ------------------------------------------------------------
static std::vector<double> phiCentersDeg() {
    std::vector<double> v(N_PHI_BINS);
    const double step = 360.0 / static_cast<double>(N_PHI_BINS);
    for (int i = 0; i < N_PHI_BINS; ++i) v[i] = (i + 0.5) * step;
    return v;
}

static void plot_group_counts(
    const std::string& group,
    const std::map<std::tuple<int,int,int,int>, HelCounts>& table,
    const std::vector<Binning>& binning_scheme,
    const std::vector<std::pair<double,double>>& xB_bins,
    const std::vector<std::pair<double,double>>& Q2_bins,
    const std::vector<std::pair<double,double>>& t_bins,
    const std::string& out_dir)
{
    static const std::vector<double> PHI = phiCentersDeg();
    std::vector<double> X(N_PHI_BINS), ex(N_PHI_BINS, 0.0);
    for (int i=0;i<N_PHI_BINS;++i) X[i] = PHI[i];
    std::error_code ec;
    std::filesystem::create_directories(out_dir, ec);

    for (int ix = 0; ix < (int)xB_bins.size(); ++ix) {
        std::set<std::pair<double,double>> q2set, tset;
        for (const auto& b : binning_scheme)
            if (std::make_pair(b.xBmin,b.xBmax)==xB_bins[ix]) {
                q2set.emplace(b.Q2min,b.Q2max);
                tset.emplace(b.tmin,b.tmax);
            }
        std::vector<std::pair<double,double>> Q2s(q2set.begin(),q2set.end());
        std::vector<std::pair<double,double>> Ts(tset.begin(),tset.end());
        if (Q2s.empty()||Ts.empty()) continue;

        const int nrows = Q2s.size();
        const int ncols = Ts.size();
        const int W=260*ncols+120, H=220*nrows+140;

        TCanvas* c = new TCanvas(Form("c_counts_%s_xB%d",group.c_str(),ix), "", W,H);
        TPad* pTop = new TPad("pTop","pTop",0.0,0.94,1.0,1.0);
        pTop->SetFillStyle(0); pTop->Draw();
        TPad* pGrid= new TPad("pGrid","pGrid",0.0,0.0,1.0,0.94);
        pGrid->SetFillStyle(0); pGrid->Draw();
        pGrid->Divide(ncols,nrows,0.0001,0.0001);

        pTop->cd();
        TLatex head; head.SetNDC(); head.SetTextSize(0.55);
        head.SetTextAlign(22);
        head.DrawLatex(0.5,0.55,Form("%s   x_{B} [%.3g,%.3g]",group.c_str(),xB_bins[ix].first,xB_bins[ix].second));

        for (int r=0;r<nrows;++r){
            for (int ccol=0;ccol<ncols;++ccol){
                pGrid->cd(r*ncols+ccol+1);
                gPad->SetGrid(1,1);
                TH1* frame=gPad->DrawFrame(0,0,360,1);
                frame->GetXaxis()->SetTitle("#phi (deg)");
                frame->GetYaxis()->SetTitle("Counts");
                std::vector<double> Yp(N_PHI_BINS,0),Ym(N_PHI_BINS,0);
                for(int ip=0;ip<N_PHI_BINS;++ip){
                    auto it=table.find({ix,r,ccol,ip});
                    if(it==table.end())continue;
                    Yp[ip]=it->second.plus;
                    Ym[ip]=it->second.minus;
                }
                TGraph* gp=new TGraph(N_PHI_BINS,X.data(),Yp.data());
                gp->SetMarkerStyle(24); gp->SetMarkerColor(kRed);
                gp->Draw("LP SAME");
                TGraph* gm=new TGraph(N_PHI_BINS,X.data(),Ym.data());
                gm->SetMarkerStyle(20); gm->SetMarkerColor(kBlue);
                gm->Draw("LP SAME");
            }
        }
        c->SaveAs((out_dir+"/plot_total_counts_"+group+"_xB_"+std::to_string(ix)+".png").c_str());
        delete c;
    }
}

// ------------------------------------------------------------
// Main compute function
// ------------------------------------------------------------
void compute_total_counts(
    const std::vector<std::string>& periods,
    const std::vector<std::string>& topologies,
    const std::vector<Binning>& binning_scheme,
    const std::map<std::string, TTree*>& dataTrees,
    const std::string& combined_cuts_json,
    const std::string& out_json_path,
    const std::string& out_root_dir)
{
    namespace fs = std::filesystem;
    const auto xB_bins = uniqueRanges(binning_scheme,'x');
    const auto Q2_bins = uniqueRanges(binning_scheme,'Q');
    const auto t_bins  = uniqueRanges(binning_scheme,'t');

    std::map<std::string, std::map<std::tuple<int,int,int,int>, HelCounts>> allGroups;

    for (const auto& period : periods) {
        auto it = dataTrees.find(toLower(period));
        if (it==dataTrees.end() || !it->second){
            std::cerr<<"[total_counts][WARN] Missing data tree for "<<period<<"\n";
            continue;
        }
        TTree* t=it->second;
        int helicity=0; double x=0,Q2=0,t1=0,phi2=0; 
        t->SetBranchAddress("helicity",&helicity);
        t->SetBranchAddress("x",&x);
        t->SetBranchAddress("Q2",&Q2);
        t->SetBranchAddress("t1",&t1);
        if(t->GetBranch("phi2")) t->SetBranchAddress("phi2",&phi2);

        std::map<std::tuple<int,int,int,int>,HelCounts> table;
        const Long64_t nent=t->GetEntries();
        for(Long64_t i=0;i<nent;++i){
            t->GetEntry(i);
            if(helicity!=+1 && helicity!=-1) continue;
            if(!std::isfinite(x)||!std::isfinite(Q2)||!std::isfinite(t1)||!std::isfinite(phi2)) continue;
            int ix=findBin(x,xB_bins), iQ=findBin(Q2,Q2_bins), itb=findBin(std::fabs(t1),t_bins), ip=phiToBin(phi2);
            if(ix<0||iQ<0||itb<0||ip<0) continue;
            auto& hc=table[{ix,iQ,itb,ip}];
            if(helicity==+1) hc.plus++; else hc.minus++;
        }

        allGroups[period]=table;

        // Plot per group
        const fs::path plot_dir = fs::path(out_root_dir)/"total_counts_plots"/period;
        plot_group_counts(period,table,binning_scheme,xB_bins,Q2_bins,t_bins,plot_dir.string());
    }

    // Write master JSON
    write_total_counts_master_json(out_json_path, allGroups, N_PHI_BINS, xB_bins, Q2_bins, t_bins);

    // Write one flat JSON per group
    const fs::path per_dir = fs::path(out_root_dir)/"jsons";
    std::error_code ec; fs::create_directories(per_dir,ec);
    for(const auto& gkv: allGroups){
        std::string fname = (per_dir/("total_counts_"+gkv.first+".json")).string();
        write_total_counts_group_json(fname,gkv.second,N_PHI_BINS,xB_bins,Q2_bins,t_bins);
    }
}