// pi0_contamination.cpp
// ------------------------------------------------------------
// Helicity-resolved π0 contamination estimator (data-driven):
//   c_h = (N_bkg_mc / N_dvcs_data_h) * (N_pi0_data_h / N_pi0_rec_mc)
// with Poisson error propagation on multiplicative factors.
//
// Inputs:
//   - DVCS DATA trees (helicity)
//   - eπ0 DATA trees (helicity)
//   - eπ0 RECO MC trees (no helicity)
//   - eπ0→DVCS BKG MC trees (no helicity)
//   - combined_cuts.json (optional 3σ exclusivity cuts on DATA/BKG hypotheses)
//
// Outputs:
//   (1) Per-period JSONs:  output/jsons/contamination/contamination_<period>.json
//   (2) Combined JSON:     output/jsons/pi0_contamination_combined.json
//       contains every period + Spring2018 + Fall2018 + 10.6_GeV
//   (3) Plots per group (periods & combined): output/contamination_plots/<group>/...
//
// Notes:
//   - No fallback logic, no copying/supplement aliasing.
//   - Topology & simple kinematic preselection applied consistently.
//   - Binning follows (xB, Q², |t|, φ) using provided binning scheme.
// ------------------------------------------------------------

#include "pi0_contamination.h"

#include <TTree.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TPad.h>
#include <TH1.h>
#include <TStyle.h>
#include <TString.h>

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

// ------------------------------------------------------------
// Config
// ------------------------------------------------------------
static constexpr int    N_PHI_BINS = 12;
static constexpr double TWO_PI     = 2.0 * M_PI;
static bool ENABLE_PLOTS = true;

// ------------------------------------------------------------
// Helpers
// ------------------------------------------------------------
static inline std::string toLower(std::string s){
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c){ return std::tolower(c); });
    return s;
}
static inline bool isSpring18(const std::string& p){
    std::string s = toLower(p);
    return s.find("sp18") != std::string::npos || s.find("spring2018") != std::string::npos;
}
static inline bool isFall18(const std::string& p){
    std::string s = toLower(p);
    return s.find("fa18") != std::string::npos || s.find("fall2018") != std::string::npos;
}

// "DVCS_Fa18_inb" → "fa18_inb"
static inline std::string periodToRunTagKey(const std::string& period){
    auto pos = period.find('_');
    if (pos==std::string::npos || pos+1>=period.size()) return toLower(period);
    return toLower(period.substr(pos+1));
}

static inline std::string topoToKey(const std::string& topoStr){
    if (topoStr == "(FD,FD)") return "FD_FD";
    if (topoStr == "(CD,FD)") return "CD_FD";
    if (topoStr == "(CD,FT)") return "CD_FT";
    return "FD_FD";
}

static inline bool passesTopology_simple(int detector1, int detector2, const std::string& topoStr){
    if (topoStr == "(FD,FD)") return (detector1==1 && detector2==1);
    if (topoStr == "(CD,FD)") return (detector1==2 && detector2==1);
    if (topoStr == "(CD,FT)") return (detector1==2 && detector2==0);
    return false;
}

static inline bool applyKinematicCuts_simple(double t1, double open_angle_ep2, double pTmiss){
    if (open_angle_ep2 <= 5.0) return false;
    if ((-t1) > 1.0)          return false;
    if (pTmiss > 0.20)        return false;
    return true;
}

static inline double wrapToTwoPi(double phi){
    if (!std::isfinite(phi)) return 0.0;
    double w = std::fmod(phi, TWO_PI);
    if (w < 0.0) w += TWO_PI;
    if (w >= TWO_PI) w = std::nextafter(TWO_PI, 0.0);
    return w;
}
static inline int phiToBin(double phi){
    double w = wrapToTwoPi(phi);
    double width = TWO_PI / static_cast<double>(N_PHI_BINS);
    int idx = static_cast<int>(std::floor(w / width));
    if (idx < 0) idx = 0;
    if (idx >= N_PHI_BINS) idx = N_PHI_BINS - 1;
    return idx;
}

static std::vector<std::pair<double,double>>
uniqueRanges(const std::vector<Binning>& scheme, char which){
    std::set<std::pair<double,double>> s;
    for (const auto& b : scheme) {
        if (which=='x') s.emplace(b.xBmin,b.xBmax);
        else if (which=='Q') s.emplace(b.Q2min,b.Q2max);
        else if (which=='t') s.emplace(b.tmin,b.tmax);
    }
    return {s.begin(), s.end()};
}

static int findBin(double v, const std::vector<std::pair<double,double>>& ranges){
    for (int i=0;i<(int)ranges.size();++i)
        if (v >= ranges[i].first && v < ranges[i].second) return i;
    return -1;
}

// ------------------------------------------------------------
// 3σ cuts loader (lightweight tolerant parser)
// ------------------------------------------------------------
struct Stats { double mean=0.0, std=0.0; };
using VarCutMap = std::map<std::string, Stats>;
using PeriodTopoCuts = std::map<std::string, VarCutMap>;  // key: "DVCS_Fa18_inb_FD_FD"

static bool loadCombinedCuts(const std::string& path, PeriodTopoCuts& out){
    std::ifstream ifs(path);
    if (!ifs) {
        std::cerr << "[pi0_contam][WARN] cannot open cuts JSON: " << path << "\n";
        return false;
    }
    std::string s((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());

    size_t pos=0;
    while (true) {
        size_t kS = s.find('"', pos); if (kS==std::string::npos) break;
        size_t kE = s.find('"', kS+1); if (kE==std::string::npos) break;
        std::string key = s.substr(kS+1, kE-kS-1);

        size_t dpos = s.find("\"data\"", kE); if (dpos==std::string::npos){ pos=kE+1; continue; }
        size_t bS = s.find('{', dpos); if (bS==std::string::npos){ pos=kE+1; continue; }

        int depth=0; size_t i=bS;
        for (; i<s.size(); ++i) { if(s[i]=='{') depth++; else if(s[i]=='}'){ depth--; if(!depth){ ++i; break; } } }
        std::string data = s.substr(bS, i-bS);

        VarCutMap cuts;
        size_t vpos=0;
        while (true) {
            size_t vS = data.find('"', vpos); if (vS==std::string::npos) break;
            size_t vE = data.find('"', vS+1); if (vE==std::string::npos) break;
            std::string var = data.substr(vS+1, vE-vS-1);

            auto readNum = [&](const char* tag)->double{
                size_t p=data.find(tag, vE); if(p==std::string::npos) return 0.0;
                p=data.find(':', p); if(p==std::string::npos) return 0.0;
                size_t a=p+1; while (a<data.size() && std::isspace((unsigned char)data[a])) ++a;
                size_t b=a;
                auto isnum=[](char c){ return std::isdigit((unsigned char)c)||c=='-'||c=='+'||c=='.'||c=='e'||c=='E'; };
                while (b<data.size() && isnum(data[b])) ++b;
                try{ return std::stod(data.substr(a,b-a)); }catch(...){ return 0.0; }
            };

            double m = readNum("\"mean\"");
            double sd= readNum("\"std\"");
            cuts[var] = Stats{m,sd};
            vpos = vE+1;
        }

        if (key.rfind("DVCS_", 0) == 0) out[key] = cuts;
        pos = kE+1;
    }
    return !out.empty();
}

static inline bool within3Sigma(double v, const Stats& s){
    return (v >= s.mean - 3.0*s.std) && (v <= s.mean + 3.0*s.std);
}
static bool passes3SigmaCuts(const VarCutMap& cuts, const std::map<std::string,double>& values){
    for (const auto& kv : cuts) {
        auto it = values.find(kv.first);
        if (it == values.end()) continue;
        if (!within3Sigma(it->second, kv.second)) return false;
    }
    return true;
}

// ------------------------------------------------------------
// Branch binders
// ------------------------------------------------------------
struct BranchBinderDVCS {
    int detector1=0, detector2=0, helicity=0;
    double t1=0, open_angle_ep2=0, pTmiss=0;
    double x=0, Q2=0, phi2=0, Delta_phi=0;

    // exclusivity-like quantities for 3σ on DVCS hypothesis
    double Emiss2=0, Mx2=0, Mx2_1=0, Mx2_2=0, theta_gamma_gamma=0, xF=0;

    bool has_d1=false, has_d2=false, has_hel=false;
    bool has_t1=false, has_oa=false, has_pT=false;
    bool has_x=false, has_Q2=false, has_phi2=false, has_Dp=false;
    bool has_Em=false, has_Mx2=false, has_Mx21=false, has_Mx22=false, has_tgg=false, has_xFv=false;

    void bind(TTree* t){
        if(!t) return;
        auto bindI=[&](const char*n,int*a,bool&f){ if(t->GetBranch(n)){ t->SetBranchAddress(n,a); f=true; } };
        auto bindD=[&](const char*n,double*a,bool&f){ if(t->GetBranch(n)){ t->SetBranchAddress(n,a); f=true; } };
        bindI("detector1",&detector1,has_d1);
        bindI("detector2",&detector2,has_d2);
        bindI("helicity",&helicity,has_hel);
        bindD("t1",&t1,has_t1);
        bindD("open_angle_ep2",&open_angle_ep2,has_oa);
        bindD("pTmiss",&pTmiss,has_pT);
        bindD("x",&x,has_x);
        bindD("Q2",&Q2,has_Q2);
        bindD("phi2",&phi2,has_phi2);
        bindD("Delta_phi",&Delta_phi,has_Dp);
        bindD("Emiss2",&Emiss2,has_Em);
        bindD("Mx2",&Mx2,has_Mx2);
        bindD("Mx2_1",&Mx2_1,has_Mx21);
        bindD("Mx2_2",&Mx2_2,has_Mx22);
        bindD("theta_gamma_gamma",&theta_gamma_gamma,has_tgg);
        bindD("xF",&xF,has_xFv);
    }
    bool readyCuts() const { return has_d1 && has_d2 && has_t1 && has_oa && has_pT && has_hel; }
    bool readyBins() const { return has_x && has_Q2 && (has_phi2 || has_Dp); }
    double phi() const { return has_phi2 ? phi2 : (has_Dp ? Delta_phi : std::numeric_limits<double>::quiet_NaN()); }
    std::map<std::string,double> cutVals() const {
        std::map<std::string,double> v;
        if (has_Dp)  v["Delta_phi"]=Delta_phi;
        if (has_tgg) v["theta_gamma_gamma"]=theta_gamma_gamma;
        if (has_pT)  v["pTmiss"]=pTmiss;
        if (has_xFv) v["xF"]=xF;
        if (has_Em)  v["Emiss2"]=Emiss2;
        if (has_Mx2) v["Mx2"]=Mx2;
        if (has_Mx21) v["Mx2_1"]=Mx2_1;
        if (has_Mx22) v["Mx2_2"]=Mx2_2;
        return v;
    }
};

struct BranchBinderEPPI0Data { // has helicity
    int detector1=0, detector2=0, helicity=0;
    double t1=0, open_angle_ep2=0, pTmiss=0;
    double x=0, Q2=0, phi2=0, Delta_phi=0;
    // exclusivity vars for eπ0 data
    double Emiss2=0, Mx2=0, Mx2_1=0, Mx2_2=0, theta_pi0_pi0=0, xF=0;

    bool has_d1=false, has_d2=false, has_hel=false;
    bool has_t1=false, has_oa=false, has_pT=false;
    bool has_x=false, has_Q2=false, has_phi2=false, has_Dp=false;
    bool has_Em=false, has_Mx2=false, has_Mx21=false, has_Mx22=false, has_tpp=false, has_xFv=false;

    void bind(TTree* t){
        if(!t) return;
        auto bindI=[&](const char*n,int*a,bool&f){ if(t->GetBranch(n)){ t->SetBranchAddress(n,a); f=true; } };
        auto bindD=[&](const char*n,double*a,bool&f){ if(t->GetBranch(n)){ t->SetBranchAddress(n,a); f=true; } };
        bindI("detector1",&detector1,has_d1);
        bindI("detector2",&detector2,has_d2);
        bindI("helicity",&helicity,has_hel);
        bindD("t1",&t1,has_t1);
        bindD("open_angle_ep2",&open_angle_ep2,has_oa);
        bindD("pTmiss",&pTmiss,has_pT);
        bindD("x",&x,has_x);
        bindD("Q2",&Q2,has_Q2);
        bindD("phi2",&phi2,has_phi2);
        bindD("Delta_phi",&Delta_phi,has_Dp);
        bindD("Emiss2",&Emiss2,has_Em);
        bindD("Mx2",&Mx2,has_Mx2);
        bindD("Mx2_1",&Mx2_1,has_Mx21);
        bindD("Mx2_2",&Mx2_2,has_Mx22);
        bindD("theta_pi0_pi0",&theta_pi0_pi0,has_tpp);
        bindD("xF",&xF,has_xFv);
    }
    bool readyCuts() const { return has_d1 && has_d2 && has_t1 && has_oa && has_pT && has_hel; }
    bool readyBins() const { return has_x && has_Q2 && (has_phi2 || has_Dp); }
    double phi() const { return has_phi2 ? phi2 : (has_Dp ? Delta_phi : std::numeric_limits<double>::quiet_NaN()); }
};

struct BranchBinderEPPI0MC { // no helicity
    int detector1=0, detector2=0;
    double t1=0, open_angle_ep2=0, pTmiss=0;
    double x=0, Q2=0, phi2=0, Delta_phi=0;
    // exclusivity vars for DVCS- and eπ0-hypothesis checks
    double Emiss2=0, Mx2=0, Mx2_1=0, Mx2_2=0, theta_pi0_pi0=0, theta_gamma_gamma=0, xF=0;

    bool has_d1=false, has_d2=false;
    bool has_t1=false, has_oa=false, has_pT=false;
    bool has_x=false, has_Q2=false, has_phi2=false, has_Dp=false;
    bool has_Em=false, has_Mx2=false, has_Mx21=false, has_Mx22=false, has_tpp=false, has_tgg=false, has_xFv=false;

    void bind(TTree* t){
        if(!t) return;
        auto bindI=[&](const char*n,int*a,bool&f){ if(t->GetBranch(n)){ t->SetBranchAddress(n,a); f=true; } };
        auto bindD=[&](const char*n,double*a,bool&f){ if(t->GetBranch(n)){ t->SetBranchAddress(n,a); f=true; } };
        bindI("detector1",&detector1,has_d1);
        bindI("detector2",&detector2,has_d2);
        bindD("t1",&t1,has_t1);
        bindD("open_angle_ep2",&open_angle_ep2,has_oa);
        bindD("pTmiss",&pTmiss,has_pT);
        bindD("x",&x,has_x);
        bindD("Q2",&Q2,has_Q2);
        bindD("phi2",&phi2,has_phi2);
        bindD("Delta_phi",&Delta_phi,has_Dp);
        bindD("Emiss2",&Emiss2,has_Em);
        bindD("Mx2",&Mx2,has_Mx2);
        bindD("Mx2_1",&Mx2_1,has_Mx21);
        bindD("Mx2_2",&Mx2_2,has_Mx22);
        bindD("theta_pi0_pi0",&theta_pi0_pi0,has_tpp);
        bindD("theta_gamma_gamma",&theta_gamma_gamma,has_tgg);
        bindD("xF",&xF,has_xFv);
    }
    bool readyCuts() const { return has_d1 && has_d2 && has_t1 && has_oa && has_pT && has_x && has_Q2 && (has_phi2 || has_Dp); }
    double phi() const { return has_phi2 ? phi2 : (has_Dp ? Delta_phi : std::numeric_limits<double>::quiet_NaN()); }

    std::map<std::string,double> cutValsForDVCS() const {
        std::map<std::string,double> v;
        if (has_Dp)  v["Delta_phi"]=Delta_phi;
        if (has_tgg) v["theta_gamma_gamma"]=theta_gamma_gamma;
        if (has_pT)  v["pTmiss"]=pTmiss;
        if (has_xFv) v["xF"]=xF;
        if (has_Em)  v["Emiss2"]=Emiss2;
        if (has_Mx2) v["Mx2"]=Mx2;
        if (has_Mx21) v["Mx2_1"]=Mx2_1;
        if (has_Mx22) v["Mx2_2"]=Mx2_2;
        return v;
    }
};

// ------------------------------------------------------------
// Containers
// ------------------------------------------------------------
using BinKey = std::tuple<int,int,int,int>; // (ix,iQ2,it,ip)
struct HelCounts { long long plus=0, minus=0; };

struct BinCounts {
    HelCounts N_data;        // DVCS data, helicity-resolved
    HelCounts N_pi0_exp;     // eπ0 data, helicity-resolved
    long long N_pi0_mc   = 0; // bkg MC (mis-ID to DVCS) -- no helicity
    long long N_pi0_reco = 0; // reco MC (true eπ0)      -- no helicity

    double c_plus = 0.0, c_plus_err = 0.0;
    double c_minus= 0.0, c_minus_err= 0.0;
};

static inline std::string keyStr(const BinKey& k){
    int a,b,c,d; std::tie(a,b,c,d)=k;
    std::ostringstream os; os<<"("<<a<<","<<b<<","<<c<<","<<d<<")";
    return os.str();
}

// Sums raw counts (not contamination) for combination
static void add_into(std::map<BinKey, BinCounts>& dst, const std::map<BinKey, BinCounts>& src){
    for (const auto& kv : src) {
        auto& d = dst[kv.first];
        d.N_data.plus     += kv.second.N_data.plus;
        d.N_data.minus    += kv.second.N_data.minus;
        d.N_pi0_exp.plus  += kv.second.N_pi0_exp.plus;
        d.N_pi0_exp.minus += kv.second.N_pi0_exp.minus;
        d.N_pi0_mc        += kv.second.N_pi0_mc;
        d.N_pi0_reco      += kv.second.N_pi0_reco;
    }
}

// Recompute c_h and σ from counts in-place
static void finalize_contamination(std::map<BinKey, BinCounts>& table){
    for (auto& kv : table) {
        auto& bc = kv.second;
        auto one=[&](long long Nexp_h, long long Ndata_h, double& c, double& e){
            if (Ndata_h<=0 || bc.N_pi0_reco<=0 || bc.N_pi0_mc<=0 || Nexp_h<=0) { c=0.0; e=0.0; return; }
            const double mc = (double)bc.N_pi0_mc;
            const double dt = (double)Ndata_h;
            const double ex = (double)Nexp_h;
            const double rc = (double)bc.N_pi0_reco;
            const double val = (mc/dt) * (ex/rc);
            auto rel=[&](double n){ return (n>0)? 1.0/std::sqrt(n) : 0.0; };
            const double rel2 = std::pow(rel(mc),2)+std::pow(rel(dt),2)+std::pow(rel(ex),2)+std::pow(rel(rc),2);
            c = val; e = val * std::sqrt(rel2);
        };
        one(bc.N_pi0_exp.plus,  bc.N_data.plus,  bc.c_plus,  bc.c_plus_err);
        one(bc.N_pi0_exp.minus, bc.N_data.minus, bc.c_minus, bc.c_minus_err);
    }
}

// ------------------------------------------------------------
// Writers
// ------------------------------------------------------------
static void write_per_period_json(
    const std::string& out_path,
    const std::map<BinKey, BinCounts>& table,
    int nPhi,
    const std::vector<std::pair<double,double>>& xB_bins,
    const std::vector<std::pair<double,double>>& Q2_bins,
    const std::vector<std::pair<double,double>>& t_bins)
{
    std::ofstream ofs(out_path);
    if (!ofs) { std::cerr<<"[pi0_contam][ERROR] open "<<out_path<<" failed\n"; return; }
    ofs<<std::fixed<<std::setprecision(8);
    ofs<<"{\n";
    ofs<<"  \"binning_meta\": {\"phi_bins\": "<<nPhi
       <<", \"xB_bins\": "<<xB_bins.size()
       <<", \"Q2_bins\": "<<Q2_bins.size()
       <<", \"t_bins\": "<<t_bins.size()<<"},\n";
    ofs<<"  \"bins\": {\n";
    bool first=true;
    for (const auto& kv : table) {
        if(!first) ofs<<",\n"; first=false;
        const auto& bc = kv.second;
        ofs<<"    \""<<keyStr(kv.first)<<"\": {"
           <<"\"N_data\":{\"helicity\":{\"+1\":"<<bc.N_data.plus<<",\"-1\":"<<bc.N_data.minus<<"},\"total\":"<<(bc.N_data.plus+bc.N_data.minus)<<"},"
           <<"\"N_pi0_exp\":{\"helicity\":{\"+1\":"<<bc.N_pi0_exp.plus<<",\"-1\":"<<bc.N_pi0_exp.minus<<"},\"total\":"<<(bc.N_pi0_exp.plus+bc.N_pi0_exp.minus)<<"},"
           <<"\"N_pi0_mc\":"<<bc.N_pi0_mc<<","
           <<"\"N_pi0_reco\":"<<bc.N_pi0_reco<<","
           <<"\"contamination\":{"
           <<"\"+1\":{\"value\":"<<bc.c_plus<<",\"err\":"<<bc.c_plus_err<<"},"
           <<"\"-1\":{\"value\":"<<bc.c_minus<<",\"err\":"<<bc.c_minus_err<<"}"
           <<"}"
           <<"}";
    }
    ofs<<"\n  }\n}\n";
    std::cout << "[pi0_contam] Wrote " << out_path << "\n";
}

static void write_combined_json(
    const std::string& out_path,
    const std::map<std::string, std::map<BinKey, BinCounts>>& groups,
    int nPhi,
    const std::vector<std::pair<double,double>>& xB_bins,
    const std::vector<std::pair<double,double>>& Q2_bins,
    const std::vector<std::pair<double,double>>& t_bins)
{
    std::ofstream ofs(out_path);
    if (!ofs) { std::cerr<<"[pi0_contam][ERROR] open "<<out_path<<" failed\n"; return; }
    ofs<<std::fixed<<std::setprecision(8);
    ofs<<"{\n";
    ofs<<"  \"binning_meta\": {\"phi_bins\": "<<nPhi
       <<", \"xB_bins\": "<<xB_bins.size()
       <<", \"Q2_bins\": "<<Q2_bins.size()
       <<", \"t_bins\": "<<t_bins.size()<<"},\n";
    ofs<<"  \"periods\": {\n";
    bool firstG=true;
    for (const auto& gkv : groups) {
        if(!firstG) ofs<<",\n"; firstG=false;
        ofs<<"    \""<<gkv.first<<"\": {\n";
        ofs<<"      \"bins\": {\n";
        bool firstB=true;
        for (const auto& kv : gkv.second) {
            if(!firstB) ofs<<",\n"; firstB=false;
            const auto& bc = kv.second;
            ofs<<"        \""<<keyStr(kv.first)<<"\": {"
               <<"\"N_data\":{\"helicity\":{\"+1\":"<<bc.N_data.plus<<",\"-1\":"<<bc.N_data.minus<<"},\"total\":"<<(bc.N_data.plus+bc.N_data.minus)<<"},"
               <<"\"N_pi0_exp\":{\"helicity\":{\"+1\":"<<bc.N_pi0_exp.plus<<",\"-1\":"<<bc.N_pi0_exp.minus<<"},\"total\":"<<(bc.N_pi0_exp.plus+bc.N_pi0_exp.minus)<<"},"
               <<"\"N_pi0_mc\":"<<bc.N_pi0_mc<<","
               <<"\"N_pi0_reco\":"<<bc.N_pi0_reco<<","
               <<"\"contamination\":{"
               <<"\"+1\":{\"value\":"<<bc.c_plus<<",\"err\":"<<bc.c_plus_err<<"},"
               <<"\"-1\":{\"value\":"<<bc.c_minus<<",\"err\":"<<bc.c_minus_err<<"}"
               <<"}"
               <<"}";
        }
        ofs<<"\n      }\n";
        ofs<<"    }";
    }
    ofs<<"\n  }\n}\n";
    std::cout << "[pi0_contam] Wrote combined " << out_path << "\n";
}

// ------------------------------------------------------------
// Plotting
// ------------------------------------------------------------
static std::vector<double> phiCentersDeg(){
    std::vector<double> v(N_PHI_BINS);
    const double step = 360.0 / double(N_PHI_BINS);
    for (int i=0;i<N_PHI_BINS;++i) v[i] = (i+0.5)*step;
    return v;
}

static void uniqueQT_for_xB(
    const std::vector<Binning>& scheme,
    const std::pair<double,double>& xBrange,
    std::vector<std::pair<double,double>>& Q2_list,
    std::vector<std::pair<double,double>>& t_list)
{
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

static int findIndex(const std::pair<double,double>& r,
                     const std::vector<std::pair<double,double>>& R){
    for (int i=0;i<(int)R.size();++i) if (R[i]==r) return i;
    return -1;
}

static void plot_group(
    const std::string& name,
    const std::map<BinKey, BinCounts>& table,
    const std::vector<Binning>& binning_scheme,
    const std::vector<std::pair<double,double>>& xB_bins,
    const std::vector<std::pair<double,double>>& Q2_bins,
    const std::vector<std::pair<double,double>>& t_bins,
    const std::string& out_dir)
{
    if (!ENABLE_PLOTS) return;
    std::error_code ec;
    std::filesystem::create_directories(out_dir, ec);

    const auto PHI = phiCentersDeg();
    std::vector<double> X(N_PHI_BINS), ex(N_PHI_BINS, 0.0);
    for (int i=0;i<N_PHI_BINS;++i) X[i]=PHI[i];

    for (int ix=0; ix<(int)xB_bins.size(); ++ix) {
        std::vector<std::pair<double,double>> Q2s, Ts;
        uniqueQT_for_xB(binning_scheme, xB_bins[ix], Q2s, Ts);
        if (Q2s.empty() || Ts.empty()) continue;

        const int nrows = (int)Q2s.size();
        const int ncols = (int)Ts.size();
        const int W = 260*ncols + 120;
        const int H = 220*nrows + 140;

        TCanvas* c = new TCanvas(Form("c_contam_%s_xB%d", name.c_str(), ix), "", W, H);

        TPad* pTop  = new TPad("pTop","pTop", 0.0, 0.94, 1.0, 1.0);
        pTop->SetFillStyle(0); pTop->SetBorderSize(0); pTop->Draw();
        TPad* pGrid = new TPad("pGrid","pGrid", 0.0, 0.0, 1.0, 0.94);
        pGrid->SetFillStyle(0); pGrid->SetBorderSize(0); pGrid->Draw();
        pGrid->cd(); pGrid->Divide(ncols, nrows, 0.0001, 0.0001);

        pTop->cd();
        TLatex head; head.SetNDC(); head.SetTextSize(0.55); head.SetTextAlign(22);
        head.DrawLatex(0.5, 0.55, Form("%s   x_{B} #in [%.3g, %.3g]", name.c_str(), xB_bins[ix].first, xB_bins[ix].second));

        for (int r=0; r<nrows; ++r) {
            int iQ = findIndex(Q2s[r], Q2_bins); if (iQ<0) continue;
            for (int cc=0; cc<ncols; ++cc) {
                int it = findIndex(Ts[cc], t_bins); if (it<0) continue;

                pGrid->cd(r*ncols + cc + 1);
                gPad->SetGrid(1,1);
                gPad->SetTopMargin(0.18);
                gPad->SetBottomMargin(0.14);
                gPad->SetLeftMargin(0.12);
                gPad->SetRightMargin(0.06);

                std::vector<double> Yp(N_PHI_BINS,0.0), Ym(N_PHI_BINS,0.0);
                std::vector<double> eYp(N_PHI_BINS,0.0), eYm(N_PHI_BINS,0.0);
                double ymax = 0.0;

                for (int ip=0; ip<N_PHI_BINS; ++ip) {
                    auto itb = table.find(BinKey(ix,iQ,it,ip));
                    if (itb == table.end()) continue;
                    const auto& bc = itb->second;
                    Yp[ip]  = bc.c_plus;
                    Ym[ip]  = bc.c_minus;
                    eYp[ip] = bc.c_plus_err;
                    eYm[ip] = bc.c_minus_err;
                    ymax = std::max(ymax, std::max(Yp[ip]+eYp[ip], Ym[ip]+eYm[ip]));
                }

                const double ymin = 0.0;
                const double ymaxPlot = std::min(1.0, (ymax>0.0 ? ymax*1.25 : 0.10));
                TH1* frame = gPad->DrawFrame(0.0, ymin, 360.0, ymaxPlot);
                frame->GetXaxis()->SetTitle("#phi (deg)");
                frame->GetYaxis()->SetTitle("#pi^{0} contamination");
                frame->GetXaxis()->CenterTitle();
                frame->GetYaxis()->CenterTitle();
                frame->GetXaxis()->SetNdivisions(505);

                TGraphErrors* grP = new TGraphErrors(N_PHI_BINS, X.data(), Yp.data(), ex.data(), eYp.data());
                TGraphErrors* grM = new TGraphErrors(N_PHI_BINS, X.data(), Ym.data(), ex.data(), eYm.data());
                grP->SetMarkerStyle(24);
                grM->SetMarkerStyle(20);
                grP->SetMarkerSize(0.9); grM->SetMarkerSize(0.9);
                grP->SetLineWidth(2);    grM->SetLineWidth(2);
                grP->Draw("P SAME");
                grM->Draw("P SAME");

                // Subplot title
                TLatex sub; sub.SetNDC(); sub.SetTextSize(0.045); sub.SetTextAlign(13);
                sub.DrawLatex(0.14, 0.96,
                    Form("Q^{2} #in [%.3g, %.3g],   -t #in [%.3g, %.3g]",
                         Q2s[r].first, Q2s[r].second, Ts[cc].first, Ts[cc].second));

                if (r==0 && cc==0) {
                    TLegend* leg = new TLegend(0.58, 0.73, 0.92, 0.92);
                    leg->SetBorderSize(1); leg->SetFillStyle(0); leg->SetTextSize(0.035);
                    leg->AddEntry(grP, "helicity +1", "p");
                    leg->AddEntry(grM, "helicity -1", "p");
                    leg->Draw();
                }
            }
        }
        std::ostringstream fout;
        fout << out_dir << "/plot_contamination_" << name << "_xB_" << ix << ".png";
        c->SaveAs(fout.str().c_str());
        delete c;
    }
}

// ------------------------------------------------------------
// Core
// ------------------------------------------------------------
void compute_pi0_contamination_helicity(
    const std::vector<std::string>& periods,
    const std::vector<std::string>& topologies,
    const std::vector<Binning>& binning_scheme,
    const std::map<std::string, TTree*>& dvcsDataTrees,   // key: runTag (e.g. "fa18_inb")
    const std::map<std::string, TTree*>& eppi0DataTrees,  // key: runTag + "_eppi0"
    const std::map<std::string, TTree*>& eppi0RecMcTrees, // key: runTag + "_rec_mc"
    const std::map<std::string, TTree*>& eppi0BkgTrees,   // key: runTag + "_bkg"
    const std::string& combined_cuts_json,
    const std::string& out_root_dir)
{
    namespace fs = std::filesystem;

    const auto xB_bins = uniqueRanges(binning_scheme,'x');
    const auto Q2_bins = uniqueRanges(binning_scheme,'Q');
    const auto t_bins  = uniqueRanges(binning_scheme,'t');
    if (xB_bins.empty() || Q2_bins.empty() || t_bins.empty()) {
        std::cerr << "[pi0_contam][ERROR] empty bin ranges\n";
        return;
    }

    // Load 3σ cuts (optional)
    PeriodTopoCuts cuts;
    loadCombinedCuts(combined_cuts_json, cuts);

    auto dvcsCutsKey = [&](const std::string& runTag, const std::string& topoKey)->std::string {
        std::string cap = runTag;
        if (!cap.empty()) cap[0] = std::toupper(static_cast<unsigned char>(cap[0]));
        for (size_t i=0;i+1<cap.size();++i) if (cap[i]=='_' && i+1<cap.size())
            cap[i+1] = std::toupper(static_cast<unsigned char>(cap[i+1]));
        return std::string("DVCS_") + cap + "_" + topoKey;
    };

    // Output paths
    const fs::path root(out_root_dir);
    const fs::path jsons_dir    = root / "jsons" / "contamination";
    const fs::path combined_out = root / "jsons" / "pi0_contamination_combined.json";
    const fs::path plots_root   = root / "contamination_plots";
    std::error_code ec;
    fs::create_directories(jsons_dir, ec);
    fs::create_directories(plots_root, ec);

    // Keep all per-period tables for combined build and plots
    std::map<std::string, std::map<BinKey, BinCounts>> groupTables;

    // ---------- Per-period loops ----------
    for (const auto& period : periods) {
        const std::string runTag = periodToRunTagKey(period);
        const std::string key_dvcs     = runTag;               // dvcs data
        const std::string key_pi0_data = runTag + "_eppi0";    // eπ0 data
        const std::string key_pi0_reco = runTag + "_rec_mc";   // eπ0 reco MC
        const std::string key_pi0_bkg  = runTag + "_bkg";      // eπ0→DVCS bkg MC

        auto itD  = dvcsDataTrees.find(key_dvcs);
        auto itED = eppi0DataTrees.find(key_pi0_data);
        auto itRM = eppi0RecMcTrees.find(key_pi0_reco);
        auto itBM = eppi0BkgTrees.find(key_pi0_bkg);

        if (itD==dvcsDataTrees.end() || !itD->second) {
            std::cerr << "[pi0_contam][WARN] Missing DVCS DATA for " << period << " (key " << key_dvcs << "). Skipping.\n";
            continue;
        }

        TTree* t_dvcs = itD->second;
        TTree* t_pi0_data = (itED!=eppi0DataTrees.end())   ? itED->second : nullptr;
        TTree* t_pi0_reco = (itRM!=eppi0RecMcTrees.end())  ? itRM->second : nullptr;
        TTree* t_pi0_bkg  = (itBM!=eppi0BkgTrees.end())    ? itBM->second : nullptr;

        if (!t_pi0_bkg || !t_pi0_reco || !t_pi0_data) {
            std::cerr << "[pi0_contam][WARN] For " << period
                      << " missing components: "
                      << "bkg="<<(t_pi0_bkg? "ok":"none")<<", "
                      << "reco="<<(t_pi0_reco? "ok":"none")<<", "
                      << "eπ0 data="<<(t_pi0_data? "ok":"none")
                      << ".\n";
        }

        std::map<BinKey, BinCounts> table;

        // DVCS data → N_data^±
        if (t_dvcs) {
            BranchBinderDVCS b; b.bind(t_dvcs);
            if (b.readyCuts() && b.readyBins()) {
                const Long64_t nent = t_dvcs->GetEntries();
                for (Long64_t i=0;i<nent;++i) {
                    t_dvcs->GetEntry(i);
                    if (!applyKinematicCuts_simple(b.t1,b.open_angle_ep2,b.pTmiss)) continue;
                    if (b.helicity!=+1 && b.helicity!=-1) continue;

                    std::string usedTopoKey;
                    for (const auto& topoStr : topologies) {
                        if (passesTopology_simple(b.detector1,b.detector2,topoStr)) {
                            usedTopoKey = topoToKey(topoStr); break;
                        }
                    }
                    if (usedTopoKey.empty()) continue;

                    // 3σ cuts if available
                    VarCutMap topoCuts;
                    auto itC = cuts.find(dvcsCutsKey(runTag, usedTopoKey));
                    if (itC != cuts.end()) topoCuts = itC->second;
                    if (!topoCuts.empty()) {
                        if (!passes3SigmaCuts(topoCuts, b.cutVals())) continue;
                    }

                    double xB=b.x, Q2=b.Q2, tt=std::fabs(b.t1), phi=b.phi();
                    if (!std::isfinite(xB)||!std::isfinite(Q2)||!std::isfinite(tt)||!std::isfinite(phi)) continue;
                    int ix=findBin(xB,xB_bins), iQ=findBin(Q2,Q2_bins), itb=findBin(tt,t_bins), ip=phiToBin(phi);
                    if (ix<0||iQ<0||itb<0||ip<0||ip>=N_PHI_BINS) continue;

                    BinKey key(ix,iQ,itb,ip);
                    if (b.helicity==+1) table[key].N_data.plus++; else table[key].N_data.minus++;
                }
            } else {
                std::cerr << "[pi0_contam][WARN] DVCS tree missing branches for " << period << "\n";
            }
        }

        // eπ0→DVCS mis-ID MC → N_pi0_mc
        if (t_pi0_bkg) {
            BranchBinderEPPI0MC b; b.bind(t_pi0_bkg);
            if (b.readyCuts()) {
                const Long64_t nent = t_pi0_bkg->GetEntries();
                for (Long64_t i=0;i<nent;++i) {
                    t_pi0_bkg->GetEntry(i);
                    if (!applyKinematicCuts_simple(b.t1,b.open_angle_ep2,b.pTmiss)) continue;

                    bool match=false; std::string usedTopoKey;
                    for (const auto& topoStr : topologies) {
                        if (passesTopology_simple(b.detector1,b.detector2,topoStr)) { match=true; usedTopoKey = topoToKey(topoStr); break; }
                    }
                    if (!match) continue;

                    // optional 3σ on DVCS hypothesis variables
                    VarCutMap topoCuts;
                    auto itC = cuts.find(dvcsCutsKey(runTag, usedTopoKey));
                    if (itC != cuts.end()) topoCuts = itC->second;
                    if (!topoCuts.empty()) {
                        if (!passes3SigmaCuts(topoCuts, b.cutValsForDVCS())) continue;
                    }

                    double xB=b.x, Q2=b.Q2, tt=std::fabs(b.t1), phi=b.phi();
                    if (!std::isfinite(xB)||!std::isfinite(Q2)||!std::isfinite(tt)||!std::isfinite(phi)) continue;
                    int ix=findBin(xB,xB_bins), iQ=findBin(Q2,Q2_bins), itb=findBin(tt,t_bins), ip=phiToBin(phi);
                    if (ix<0||iQ<0||itb<0||ip<0||ip>=N_PHI_BINS) continue;

                    table[BinKey(ix,iQ,itb,ip)].N_pi0_mc++;
                }
            } else {
                std::cerr << "[pi0_contam][WARN] eπ0_bkg MC missing branches for " << period << "\n";
            }
        }

        // eπ0 reco MC → N_pi0_reco
        if (t_pi0_reco) {
            BranchBinderEPPI0MC b; b.bind(t_pi0_reco);
            if (b.readyCuts()) {
                const Long64_t nent = t_pi0_reco->GetEntries();
                for (Long64_t i=0;i<nent;++i) {
                    t_pi0_reco->GetEntry(i);
                    if (!applyKinematicCuts_simple(b.t1,b.open_angle_ep2,b.pTmiss)) continue;

                    bool match=false;
                    for (const auto& topoStr : topologies) {
                        if (passesTopology_simple(b.detector1,b.detector2,topoStr)) { match=true; break; }
                    }
                    if (!match) continue;

                    double xB=b.x, Q2=b.Q2, tt=std::fabs(b.t1), phi=b.phi();
                    if (!std::isfinite(xB)||!std::isfinite(Q2)||!std::isfinite(tt)||!std::isfinite(phi)) continue;
                    int ix=findBin(xB,xB_bins), iQ=findBin(Q2,Q2_bins), itb=findBin(tt,t_bins), ip=phiToBin(phi);
                    if (ix<0||iQ<0||itb<0||ip<0||ip>=N_PHI_BINS) continue;

                    table[BinKey(ix,iQ,itb,ip)].N_pi0_reco++;
                }
            } else {
                std::cerr << "[pi0_contam][WARN] eπ0 reco MC missing branches for " << period << "\n";
            }
        }

        // eπ0 DATA → N_pi0_exp^±
        if (t_pi0_data) {
            BranchBinderEPPI0Data b; b.bind(t_pi0_data);
            if (b.readyCuts() && b.readyBins()) {
                const Long64_t nent = t_pi0_data->GetEntries();
                for (Long64_t i=0;i<nent;++i) {
                    t_pi0_data->GetEntry(i);
                    if (!applyKinematicCuts_simple(b.t1,b.open_angle_ep2,b.pTmiss)) continue;
                    if (b.helicity!=+1 && b.helicity!=-1) continue;

                    bool match=false;
                    for (const auto& topoStr : topologies) {
                        if (passesTopology_simple(b.detector1,b.detector2,topoStr)) { match=true; break; }
                    }
                    if (!match) continue;

                    double xB=b.x, Q2=b.Q2, tt=std::fabs(b.t1), phi=b.phi();
                    if (!std::isfinite(xB)||!std::isfinite(Q2)||!std::isfinite(tt)||!std::isfinite(phi)) continue;
                    int ix=findBin(xB,xB_bins), iQ=findBin(Q2,Q2_bins), itb=findBin(tt,t_bins), ip=phiToBin(phi);
                    if (ix<0||iQ<0||itb<0||ip<0||ip>=N_PHI_BINS) continue;

                    BinKey key(ix,iQ,itb,ip);
                    if (b.helicity==+1) table[key].N_pi0_exp.plus++; else table[key].N_pi0_exp.minus++;
                }
            } else {
                std::cerr << "[pi0_contam][WARN] eπ0 DATA missing branches for " << period << "\n";
            }
        }

        // Compute contamination for this period
        finalize_contamination(table);

        // Write per-period JSON
        const std::string out_per =
            (jsons_dir / ("contamination_" + period + ".json")).string();
        write_per_period_json(out_per, table, N_PHI_BINS, xB_bins, Q2_bins, t_bins);

        // Plot per-period
        const std::string plot_dir = (plots_root / period).string();
        plot_group(period, table, binning_scheme, xB_bins, Q2_bins, t_bins, plot_dir);

        // Stash for combined
        groupTables[period] = std::move(table);
    }

    // ---------- Build combined groups from raw counts ----------
    std::map<BinKey, BinCounts> spring_sum, fall_sum, all_sum;
    for (const auto& g : groupTables) {
        const std::string& name = g.first;
        if (isSpring18(name)) add_into(spring_sum, g.second);
        if (isFall18(name))   add_into(fall_sum,   g.second);
        add_into(all_sum, g.second);
    }
    if (!spring_sum.empty()) { finalize_contamination(spring_sum); groupTables["Spring2018"] = spring_sum; }
    if (!fall_sum.empty())   { finalize_contamination(fall_sum);   groupTables["Fall2018"]   = fall_sum;   }
    if (!all_sum.empty())    { finalize_contamination(all_sum);    groupTables["10.6_GeV"]   = all_sum;    }

    // Write combined JSON with everything (periods + combined)
    write_combined_json(combined_out.string(), groupTables, N_PHI_BINS, xB_bins, Q2_bins, t_bins);

    // Plots for combined groups too
    for (const auto& g : {"Spring2018","Fall2018","10.6_GeV"}) {
        auto it = groupTables.find(g);
        if (it == groupTables.end()) continue;
        const std::string dir = (plots_root / std::string(g)).string();
        plot_group(g, it->second, binning_scheme, xB_bins, Q2_bins, t_bins, dir);
    }
}