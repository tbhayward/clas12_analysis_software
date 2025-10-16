// pi0_contamination.cpp
//
// Helicity-resolved π0 contamination (per bin & φ), with Fa18_inb_supp
// reusing the Fa18_inb MC. Writes per-period JSONs under
//   output/jsons/contamination/contamination_<period>.json
// a combined JSON:
//   output/jsons/pi0_contamination_combined.json
// and plots under:
//   output/contamination_plots/<runTag>/plot_contam_<period>_xB_<ix>.png
//
// Contamination estimate (simple and stable):
//   For each (ix,iQ,it,ip, helicity h),
//     c_h = clamp( N_bkg_mc_h / max(N_sig_mc_h, 1), 0, 0.95 )
// where "sig_mc" is the reconstructed eπ0 MC (signal) and "bkg_mc"
// is the DVCS-background MC provided for eπ0.
//
// NOTE: This file uses the same simple kinematic & topology cuts
// helpers as other stages, so everything stays consistent.

#include "pi0_contamination.h"

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TLatex.h>
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

constexpr int    N_PHI_BINS = 12;
constexpr double TWO_PI     = 2.0 * M_PI;

struct Stats { double mean=0.0; double std=0.0; };
using BinKey = std::tuple<int,int,int,int>; // (ix,iQ,it,ip)

// -------------------- style --------------------
struct StyleInit {
    StyleInit() {
        gStyle->SetOptTitle(0);
        gStyle->SetOptStat(0);
        gStyle->SetFrameLineWidth(2);
        gStyle->SetLineWidth(2);
        gStyle->SetPadTickX(1);
        gStyle->SetPadTickY(1);
        gStyle->SetLegendBorderSize(1);
        gStyle->SetTitleFont(42, "XYZ");
        gStyle->SetLabelFont(42, "XYZ");
        gStyle->SetTextFont(42);
    }
} _style_bootstrap;

// -------------------- utils --------------------
static inline std::string toLower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c){ return std::tolower(c); });
    return s;
}

// "DVCS_Fa18_inb" -> "fa18_inb"
static std::string periodToRunTagKey(const std::string& period) {
    auto pos = period.find('_');
    if (pos == std::string::npos || pos + 1 >= period.size()) return toLower(period);
    return toLower(period.substr(pos + 1));
}

static inline std::vector<std::pair<double,double>> uniqueRanges(const std::vector<Binning>& scheme, char which) {
    std::set<std::pair<double,double>> s;
    for (const auto& b : scheme) {
        if (which == 'x') s.emplace(b.xBmin, b.xBmax);
        else if (which == 'Q') s.emplace(b.Q2min, b.Q2max);
        else if (which == 't') s.emplace(b.tmin, b.tmax);
    }
    return std::vector<std::pair<double,double>>(s.begin(), s.end());
}

static inline int findBin(double v, const std::vector<std::pair<double,double>>& ranges) {
    for (int i=0;i<(int)ranges.size();++i)
        if (v>=ranges[i].first && v<ranges[i].second) return i;
    return -1;
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
    double width = TWO_PI / double(N_PHI_BINS);
    int idx = std::max(0, std::min(int(std::floor(w/width)), N_PHI_BINS-1));
    return idx;
}

static inline bool passesTopology_simple(int detector1, int detector2, const std::vector<std::string>& tops){
    for (const auto& t : tops){
        if (t=="(FD,FD)" && detector1==1 && detector2==1) return true;
        if (t=="(CD,FD)" && detector1==2 && detector2==1) return true;
        if (t=="(CD,FT)" && detector1==2 && detector2==0) return true;
    }
    return false;
}

static inline bool applyKinematicCuts_simple(double t1, double open_angle_ep2, double pTmiss) {
    if (open_angle_ep2 <= 5.0) return false;
    if ((-t1) > 1.0)          return false;
    if (pTmiss > 0.20)        return false;
    return true;
}

// For a given xB slice, list only the Q² and |t| ranges present in CSV.
static void uniqueQT_for_xB(
    const std::vector<Binning>& scheme,
    const std::pair<double,double>& xBrange,
    std::vector<std::pair<double,double>>& Q2_list,
    std::vector<std::pair<double,double>>& t_list
) {
    std::set<std::pair<double,double>> qs, ts;
    for (const auto& b : scheme) {
        if (std::make_pair(b.xBmin,b.xBmax) == xBrange) {
            qs.emplace(b.Q2min,b.Q2max);
            ts.emplace(b.tmin,b.tmax);
        }
    }
    Q2_list.assign(qs.begin(), qs.end());
    t_list.assign(ts.begin(), ts.end());
}

static inline std::vector<double> phiCentersDeg() {
    std::vector<double> d(N_PHI_BINS);
    const double step = 360.0/double(N_PHI_BINS);
    for (int i=0;i<N_PHI_BINS;++i) d[i] = (i+0.5)*step;
    return d;
}

// ----------- branches (minimal set we need) -----------
struct Branch {
    int helicity=0, detector1=0, detector2=0;
    double t1=0, open_angle_ep2=0, pTmiss=0, x=0, Q2=0, phi2=0, Delta_phi=0;
    bool hasH=false, hasTopo=false, hasCuts=false, hasBins=false, hasPhi2=false, hasDp=false;
    void bind(TTree* t){
        if (!t) return;
        auto bindI=[&](const char*n,int*a,bool&f){ if(t->GetBranch(n)){ t->SetBranchAddress(n,a); f=true; } };
        auto bindD=[&](const char*n,double*a,bool&f){ if(t->GetBranch(n)){ t->SetBranchAddress(n,a); f=true; } };
        bindI("helicity",&helicity,hasH);
        bindI("detector1",&detector1,hasTopo);
        bindI("detector2",&detector2,hasTopo);
        bindD("t1",&t1,hasCuts);
        bindD("open_angle_ep2",&open_angle_ep2,hasCuts);
        bindD("pTmiss",&pTmiss,hasCuts);
        bindD("x",&x,hasBins);
        bindD("Q2",&Q2,hasBins);
        bindD("phi2",&phi2,hasPhi2);
        bindD("Delta_phi",&Delta_phi,hasDp);
    }
    inline bool ok() const { return hasH && hasTopo && hasCuts && hasBins && (hasPhi2||hasDp); }
    inline double phi() const { return hasPhi2 ? phi2 : (hasDp ? Delta_phi : std::numeric_limits<double>::quiet_NaN()); }
};

// --- 3σ cuts: use same "combined_cuts.json" light parser as elsewhere ---
using VarCutMap = std::map<std::string, Stats>;
using PeriodTopoCuts = std::map<std::string, VarCutMap>;

static bool loadCombinedCuts(const std::string& path, PeriodTopoCuts& out) {
    std::ifstream ifs(path);
    if (!ifs) return false;
    std::string s((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    size_t pos = 0;
    while (true) {
        size_t kS = s.find('"', pos);
        if (kS == std::string::npos) break;
        size_t kE = s.find('"', kS+1);
        if (kE == std::string::npos) break;
        std::string key = s.substr(kS+1, kE-kS-1);

        size_t dpos = s.find("\"data\"", kE);
        if (dpos == std::string::npos) { pos = kE+1; continue; }
        size_t bS = s.find('{', dpos);
        if (bS == std::string::npos) { pos = kE+1; continue; }
        int depth=0; size_t i=bS;
        for (; i<s.size(); ++i) { if(s[i]=='{') depth++; else if(s[i]=='}'){ depth--; if(!depth){ ++i; break; } } }
        std::string data = s.substr(bS, i-bS);

        VarCutMap cuts;
        size_t vpos = 0;
        while (true) {
            size_t vS = data.find('"', vpos);
            if (vS == std::string::npos) break;
            size_t vE = data.find('"', vS+1);
            if (vE == std::string::npos) break;
            std::string var = data.substr(vS+1, vE-vS-1);

            auto readNum = [&](const char* tag)->double{
                size_t p = data.find(tag, vE);
                if (p == std::string::npos) return 0.0;
                p = data.find(':', p); if(p==std::string::npos) return 0.0;
                size_t j = p+1; while (j<data.size() && std::isspace((unsigned char)data[j])) ++j;
                size_t k=j;
                auto isnum=[&](char c){return (std::isdigit((unsigned char)c)||c=='-'||c=='+'||c=='.'||c=='e'||c=='E');};
                while (k<data.size() && isnum(data[k])) ++k;
                try { return std::stod(data.substr(j,k-j)); } catch (...) { return 0.0; }
            };
            double mean = readNum("\"mean\"");
            double std  = readNum("\"std\"");
            cuts[var] = Stats{mean, std};
            vpos = vE+1;
        }

        if (key.rfind("DVCS_", 0) == 0) out[key] = cuts;
        pos = kE+1;
    }
    return !out.empty();
}

static inline bool within3Sigma(double v, const Stats& s) {
    return (v >= s.mean - 3.0*s.std) && (v <= s.mean + 3.0*s.std);
}

// We only apply 3σ cuts for variables present in the tree & in the cuts map.
static bool passes3SigmaCuts(const VarCutMap& /*cuts*/, const Branch& /*b*/) {
    // Light placeholder: if you later want to map var names -> branch values exactly
    // (e.g., Delta_phi, pTmiss, theta_gamma_gamma), wire it here.
    // For now we skip applying them at the contamination stage; exclusivity already applied upstream.
    return true;
}

// ---------- MC alias for Fa18_inb_supp ----------
static inline std::string period_alias_for_mc(const std::string& period) {
    if (period == "DVCS_Fa18_inb_supp") return "DVCS_Fa18_inb";
    return period;
}

// ---------- minimal JSON writer for contamination ----------
struct ContamBin { double c_plus=0, c_plus_err=0, c_minus=0, c_minus_err=0; };

static inline std::string key4(int ix,int iQ,int it,int ip){
    std::ostringstream os; os<<"("<<ix<<","<<iQ<<","<<it<<","<<ip<<")"; return os.str();
}

static void write_period_contam_json(
    const std::string& out_path,
    int nPhi,
    const std::vector<std::pair<double,double>>& xB_bins,
    const std::vector<std::pair<double,double>>& Q2_bins,
    const std::vector<std::pair<double,double>>& t_bins,
    const std::map<BinKey, ContamBin>& cmap)
{
    std::ofstream ofs(out_path);
    if (!ofs) { std::cerr<<"[pi0_contam][ERROR] Cannot open "<<out_path<<"\n"; return; }
    ofs<<std::fixed<<std::setprecision(8);
    ofs<<"{\n";
    ofs<<"  \"binning_meta\": {\"phi_bins\": "<<nPhi<<", \"xB_bins\": "<<xB_bins.size()<<", \"Q2_bins\": "<<Q2_bins.size()<<", \"t_bins\": "<<t_bins.size()<<"},\n";
    ofs<<"  \"bins\": {\n";
    bool first=true;
    for (const auto& kv : cmap){
        if (!first) ofs<<",\n"; first=false;
        int ix,iQ,it,ip; std::tie(ix,iQ,it,ip)=kv.first;
        const auto& cb = kv.second;
        ofs<<"    \""<<key4(ix,iQ,it,ip)<<"\": {";
        ofs<<"\"contamination\":{\"+1\":{\"value\":"<<cb.c_plus<<",\"err\":"<<cb.c_plus_err<<"},"
              "\"-1\":{\"value\":"<<cb.c_minus<<",\"err\":"<<cb.c_minus_err<<"}}";
        ofs<<"}";
    }
    ofs<<"\n  }\n}\n";
}

// ---------- plotting ----------
static inline int findIndex(const std::pair<double,double>& r,
                            const std::vector<std::pair<double,double>>& R){
    for (int i=0;i<(int)R.size();++i) if (R[i]==r) return i;
    return -1;
}

static void plot_contamination_for_period(
    const std::string& period,
    const std::vector<Binning>& binning_scheme,
    const std::vector<std::pair<double,double>>& xB_bins,
    const std::vector<std::pair<double,double>>& Q2_bins,
    const std::vector<std::pair<double,double>>& t_bins,
    const std::map<BinKey, ContamBin>& cmap,
    const std::string& out_dir_plots)
{
    namespace fs = std::filesystem;
    std::error_code ec;
    fs::create_directories(out_dir_plots, ec);

    static const auto PHI_DEG = phiCentersDeg();

    for (int ix=0; ix<(int)xB_bins.size(); ++ix){
        const auto xb = xB_bins[ix];

        std::vector<std::pair<double,double>> Q2_slice, t_slice;
        uniqueQT_for_xB(binning_scheme, xb, Q2_slice, t_slice);
        if (Q2_slice.empty() || t_slice.empty()) continue;

        const int nrows = (int)t_slice.size();
        const int ncols = (int)Q2_slice.size();

        const int W = 280*ncols + 160;
        const int H = 240*nrows + 170;

        std::ostringstream cname; cname<<"c_contam_"<<period<<"_xB"<<ix;
        TCanvas* c = new TCanvas(cname.str().c_str(), cname.str().c_str(), W, H);

        // title pad
        TPad* pTop  = new TPad("pTop","pTop", 0.0, 0.915, 1.0, 1.0);
        pTop->SetFillStyle(0); pTop->SetBorderSize(0); pTop->Draw();

        // grid pad
        TPad* pGrid = new TPad("pGrid","pGrid", 0.0, 0.00, 1.0, 0.915);
        pGrid->SetFillStyle(0); pGrid->SetBorderSize(0); pGrid->Draw();
        pGrid->cd();
        pGrid->Divide(ncols, nrows, 0.0001, 0.0001);

        // title text
        pTop->cd();
        TLatex head;
        head.SetNDC(); head.SetTextAlign(22); head.SetTextFont(42); head.SetTextSize(0.36);
        head.DrawLatex(0.5, 0.55, Form("#pi^{0} contamination  %s   x_{B} #in [%.2g, %.2g]",
                                       period.c_str(), xb.first, xb.second));

        for (int r=0; r<nrows; ++r){
            const int it_glob = findIndex(t_slice[r], t_bins);
            if (it_glob<0) continue;
            for (int cc=0; cc<ncols; ++cc){
                const int iQ_glob = findIndex(Q2_slice[cc], Q2_bins);
                if (iQ_glob<0) continue;

                pGrid->cd(r*ncols + cc + 1);
                gPad->SetGrid(1,1);
                gPad->SetTopMargin(0.08);
                gPad->SetBottomMargin(0.18);
                gPad->SetLeftMargin(0.15);
                gPad->SetRightMargin(0.10);

                TH1* frame = gPad->DrawFrame(0.0, 0.0, 360.0, 0.6);
                frame->GetXaxis()->SetTitle("#phi (deg)");
                frame->GetYaxis()->SetTitle("contamination");
                frame->GetXaxis()->CenterTitle();
                frame->GetYaxis()->CenterTitle();
                frame->GetXaxis()->SetTitleSize(0.060);
                frame->GetYaxis()->SetTitleSize(0.060);
                frame->GetYaxis()->SetLabelSize(0.048);
                frame->GetXaxis()->SetTitleOffset(1.25);
                frame->GetYaxis()->SetTitleOffset(1.35);

                // build series for +1 and -1
                std::vector<double> x, yP, yM;
                for (int ip=0; ip<N_PHI_BINS; ++ip){
                    auto it = cmap.find(std::make_tuple(ix,iQ_glob,it_glob,ip));
                    if (it==cmap.end()) continue;
                    x.push_back(PHI_DEG[ip]);
                    yP.push_back(std::max(0.0, std::min(0.95, it->second.c_plus)));
                    yM.push_back(std::max(0.0, std::min(0.95, it->second.c_minus)));
                }

                if (!x.empty()){
                    TGraph* gP = new TGraph((int)x.size(), x.data(), yP.data());
                    TGraph* gM = new TGraph((int)x.size(), x.data(), yM.data());
                    gP->SetLineWidth(2); gP->SetMarkerStyle(20); gP->Draw("LP SAME");
                    gM->SetLineWidth(2); gM->SetMarkerStyle(24); gM->Draw("LP SAME");
                    TLegend* leg = new TLegend(0.50, 0.74, 0.90, 0.92);
                    leg->SetTextFont(42); leg->SetTextSize(0.040);
                    leg->AddEntry(gP, "+1 helicity", "lp");
                    leg->AddEntry(gM, "-1 helicity", "lp");
                    leg->Draw();
                }

                TLatex lab; lab.SetNDC(); lab.SetTextSize(0.040); lab.SetTextAlign(11); lab.SetTextFont(42);
                lab.DrawLatex(0.15, 0.94,
                    Form("Q^{2} #in [%.2g, %.2g],   -t #in [%.2g, %.2g]",
                         Q2_slice[cc].first, Q2_slice[cc].second,
                         t_slice[r].first,   t_slice[r].second));
            }
        }

        std::ostringstream fout;
        fout<<out_dir_plots<<"/plot_contam_"<<period<<"_xB_"<<ix<<".png";
        c->SaveAs(fout.str().c_str());
        delete c;
    }
}

} // namespace

// ======================================================================
// Public driver
// ======================================================================
void compute_pi0_contamination_helicity(
    const std::vector<std::string>& periods,
    const std::vector<std::string>& topologies,
    const std::vector<Binning>& binning_scheme,
    const std::map<std::string, TTree*>& /*dvcsDataTrees*/, // not used directly here
    const std::map<std::string, TTree*>& eppi0DataTrees,
    const std::map<std::string, TTree*>& eppi0RecMcTrees,
    const std::map<std::string, TTree*>& eppi0BkgTrees,
    const std::string& combined_cuts_json_path,
    const std::string& out_root_dir)
{
    namespace fs = std::filesystem;

    const auto xB_bins = uniqueRanges(binning_scheme, 'x');
    const auto Q2_bins = uniqueRanges(binning_scheme, 'Q');
    const auto t_bins  = uniqueRanges(binning_scheme, 't');

    PeriodTopoCuts cuts;
    if (!loadCombinedCuts(combined_cuts_json_path, cuts)) {
        std::cerr << "[pi0_contam][WARN] Could not load combined cuts: " << combined_cuts_json_path << "\n";
    }

    // Storage for combined JSON
    std::map<std::string, std::map<BinKey, ContamBin>> allPeriod;

    // Per-period loop
    for (const auto& period : periods) {

        const std::string mc_period = period_alias_for_mc(period);
        const std::string runTag_mc = periodToRunTagKey(mc_period); // fa18_inb for supp
        const std::string runTag    = periodToRunTagKey(period);    // fa18_inb_supp etc

        // Build expected keys
        const std::string reco_key = runTag_mc + "_rec_mc"; // e.g. fa18_inb_rec_mc
        const std::string bkg_key  = runTag_mc + "_bkg";    // e.g. fa18_inb_bkg

        TTree* reco_mc = nullptr;
        TTree* bkg_mc  = nullptr;

        auto itR = eppi0RecMcTrees.find(reco_key);
        if (itR != eppi0RecMcTrees.end()) reco_mc = itR->second;

        auto itB = eppi0BkgTrees.find(bkg_key);
        if (itB != eppi0BkgTrees.end()) bkg_mc = itB->second;

        if (!reco_mc || !bkg_mc) {
            std::cerr << "[pi0_contam][WARN] Missing π0 MC for '"<<period<<"'; "
                      << "bkg="<<(bkg_mc?"ok":"none")<<" reco="<<(reco_mc?"ok":"none")<<". Continuing with what is available.\n";
        }

        // Count MC per (ix,iQ,it,ip, helicity)
        std::map<BinKey, long long> mc_sig_plus, mc_sig_minus;
        std::map<BinKey, long long> mc_bkg_plus, mc_bkg_minus;

        auto fillFromTree = [&](TTree* t, std::map<BinKey,long long>& Hplus, std::map<BinKey,long long>& Hminus){
            if (!t) return;
            Branch b; b.bind(t);
            if (!b.ok()) return;
            const Long64_t nent = t->GetEntries();
            for (Long64_t i=0;i<nent;++i){
                t->GetEntry(i);
                if (!applyKinematicCuts_simple(b.t1,b.open_angle_ep2,b.pTmiss)) continue;
                if (!passesTopology_simple(b.detector1,b.detector2,topologies)) continue;
                if (b.helicity!=+1 && b.helicity!=-1) continue;

                double xB=b.x, Q2=b.Q2, tt=fabs(b.t1), phi=b.phi();
                if (!std::isfinite(xB)||!std::isfinite(Q2)||!std::isfinite(tt)||!std::isfinite(phi)) continue;
                int ix=findBin(xB,xB_bins), iQ=findBin(Q2,Q2_bins), it=findBin(tt,t_bins), ip=phiToBin(phi);
                if (ix<0||iQ<0||it<0||ip<0||ip>=N_PHI_BINS) continue;

                BinKey key(ix,iQ,it,ip);
                if (b.helicity==+1) Hplus[key]++; else Hminus[key]++;
            }
        };

        fillFromTree(reco_mc, mc_sig_plus, mc_sig_minus);
        fillFromTree(bkg_mc,  mc_bkg_plus, mc_bkg_minus);

        // Build contamination map
        std::map<BinKey, ContamBin> cmap;
        auto clamp01 = [](double v){ return std::max(0.0, std::min(0.95, v)); };

        // unify keys present in either
        std::set<BinKey> keys;
        for (const auto& kv : mc_sig_plus) keys.insert(kv.first);
        for (const auto& kv : mc_sig_minus) keys.insert(kv.first);
        for (const auto& kv : mc_bkg_plus) keys.insert(kv.first);
        for (const auto& kv : mc_bkg_minus) keys.insert(kv.first);

        for (const auto& k : keys){
            const long long sP = (mc_sig_plus.count(k)? mc_sig_plus.at(k) : 0LL);
            const long long sM = (mc_sig_minus.count(k)? mc_sig_minus.at(k) : 0LL);
            const long long bP = (mc_bkg_plus.count(k)? mc_bkg_plus.at(k) : 0LL);
            const long long bM = (mc_bkg_minus.count(k)? mc_bkg_minus.at(k) : 0LL);

            // fraction with simple binomial error approximation
            ContamBin cb;
            if (sP>0) {
                cb.c_plus = clamp01( double(bP)/double(std::max<long long>(sP,1)) );
                cb.c_plus_err = std::sqrt( cb.c_plus * (1.0 - cb.c_plus) / double(std::max<long long>(sP,1)) );
            } else { cb.c_plus = 0.0; cb.c_plus_err = 0.0; }
            if (sM>0) {
                cb.c_minus = clamp01( double(bM)/double(std::max<long long>(sM,1)) );
                cb.c_minus_err = std::sqrt( cb.c_minus * (1.0 - cb.c_minus) / double(std::max<long long>(sM,1)) );
            } else { cb.c_minus = 0.0; cb.c_minus_err = 0.0; }

            cmap[k]=cb;
        }

        // write per-period JSON
        const fs::path out_json_dir = fs::path(out_root_dir)/"jsons"/"contamination";
        std::error_code ec; fs::create_directories(out_json_dir, ec);
        const fs::path out_json = out_json_dir/("contamination_"+period+".json");
        write_period_contam_json(out_json.string(), N_PHI_BINS, xB_bins, Q2_bins, t_bins, cmap);
        std::cout<<"[pi0_contam] Wrote "<<out_json.string()<<"\n";

        // plots
        const fs::path plots_dir = fs::path(out_root_dir)/"contamination_plots"/periodToRunTagKey(period);
        plot_contamination_for_period(period, binning_scheme, xB_bins, Q2_bins, t_bins, cmap, plots_dir.string());

        allPeriod[period] = std::move(cmap);
    }

    // combined JSON (flat per-period object)
    {
        const fs::path out_comb = fs::path(out_root_dir)/"jsons"/"pi0_contamination_combined.json";
        std::ofstream ofs(out_comb.string());
        if (ofs){
            ofs<<std::fixed<<std::setprecision(8);
            ofs<<"{\n  \"periods\": {\n";
            bool firstP=true;
            for (const auto& pkv : allPeriod){
                if (!firstP) ofs<<",\n"; firstP=false;
                ofs<<"    \""<<pkv.first<<"\": {\n      \"bins\": {\n";
                bool firstB=true;
                for (const auto& kv : pkv.second){
                    if (!firstB) ofs<<",\n"; firstB=false;
                    int ix,iQ,it,ip; std::tie(ix,iQ,it,ip)=kv.first;
                    const auto& cb=kv.second;
                    ofs<<"        \""<<key4(ix,iQ,it,ip)<<"\": {"
                       <<"\"contamination\":{\"+1\":{\"value\":"<<cb.c_plus<<",\"err\":"<<cb.c_plus_err<<"},"
                       <<"\"-1\":{\"value\":"<<cb.c_minus<<",\"err\":"<<cb.c_minus_err<<"}}}";
                }
                ofs<<"\n      }\n    }";
            }
            ofs<<"\n  }\n}\n";
            std::cout<<"[pi0_contam] Wrote combined "<<out_comb.string()<<"\n";
        }
    }
}