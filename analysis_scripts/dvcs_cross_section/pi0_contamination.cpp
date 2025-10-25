// pi0_contamination.cpp
// ------------------------------------------------------------
// Helicity-resolved π0 contamination estimator (STRICT, fail-fast):
//   c_h = (N_bkg_mc / N_dvcs_data_h) * (N_pi0_data_h / N_pi0_rec_mc)
// with Poisson error propagation.
//
// Inputs (all required per period), using canonical keys that match load_trees.cpp:
//   - DVCS DATA trees (helicity)            key: DVCS_<Period>_<beam>   e.g. DVCS_Sp18_inb
//   - eπ0 DATA trees (helicity)             key: DVCS_<...>_eppi0       e.g. DVCS_Sp18_inb_eppi0
//   - eπ0 RECO MC trees (no helicity)       key: DVCS_<...>_rec_mc      e.g. DVCS_Sp18_inb_rec_mc
//   - eπ0→DVCS BKG MC trees (no helicity)   key: DVCS_<...>_bkg         e.g. DVCS_Sp18_inb_bkg
//   - combined_cuts.json   (required; 3σ exclusivity cuts on DVCS hypothesis)
//       * If a cut variable exists in JSON but branch missing in the tree -> FATAL
//
// Outputs:
//   (1) Per-period JSONs:  <out_root_dir>/jsons/contamination/contamination_<period>.json
//   (2) Combined JSON:     <out_root_dir>/jsons/pi0_contamination_combined.json
//   (3) Plots per group:   <out_root_dir>/contamination_plots/<group>/...
//
// Notes:
//   - Uses only phi2 (no Delta_phi fallback).
//   - Binning is (xB, Q², |t|, φ) from provided binning scheme.
//   - Any missing lookup / I/O failure is FATAL.
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
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace {

// ------------------------------------------------------------
// Config
// ------------------------------------------------------------
static constexpr int    N_PHI_BINS = 12;
static constexpr double TWO_PI     = 2.0 * M_PI;
static bool ENABLE_PLOTS = true;

// ------------------------------------------------------------
// Fatal helper
// ------------------------------------------------------------
[[noreturn]] static void fatal(const std::string& msg) {
    std::cerr << "[pi0_contam][FATAL] " << msg << std::endl;
    std::exit(EXIT_FAILURE);
}

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

static inline std::string topoToKey(const std::string& topoStr){
    if (topoStr == "(FD,FD)") return "FD_FD";
    if (topoStr == "(CD,FD)") return "CD_FD";
    if (topoStr == "(CD,FT)") return "CD_FT";
    fatal(std::string("Unknown topology string: ") + topoStr);
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
    if (!std::isfinite(phi)) return std::numeric_limits<double>::quiet_NaN();
    double w = std::fmod(phi, TWO_PI);
    if (w < 0.0) w += TWO_PI;
    if (w >= TWO_PI) w = std::nextafter(TWO_PI, 0.0);
    return w;
}
static inline int phiToBin(double phi){
    double w = wrapToTwoPi(phi);
    if (!std::isfinite(w)) return -1;
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
// 3σ cuts loader (STRICT)
// ------------------------------------------------------------
struct Stats { double mean=0.0, std=0.0; };
using VarCutMap = std::map<std::string, Stats>;
using PeriodTopoCuts = std::map<std::string, VarCutMap>;  // key: "DVCS_<RunTag>_<TopoKey>" (RunTag EXACT, e.g. DVCS_Sp18_inb)

static bool parseNumber(const std::string& s, size_t posColon, double& out){
    size_t a = posColon + 1;
    while (a < s.size() && std::isspace(static_cast<unsigned char>(s[a]))) ++a;
    size_t b = a;
    auto isnum=[](char c){ return std::isdigit((unsigned char)c)||c=='-'||c=='+'||c=='.'||c=='e'||c=='E'; };
    while (b < s.size() && isnum(s[b])) ++b;
    try { out = std::stod(s.substr(a, b - a)); return true; } catch (...) { return false; }
}

static void loadCombinedCuts_STRICT(const std::string& path, PeriodTopoCuts& out){
    std::ifstream ifs(path);
    if (!ifs) fatal(std::string("Cannot open cuts JSON: ") + path);
    std::string s((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    if (s.empty()) fatal(std::string("Cuts JSON is empty: ") + path);

    auto parse_num_after_colon = [](const std::string& src, size_t colonPos)->double{
        size_t a = colonPos + 1;
        while (a < src.size() && std::isspace((unsigned char)src[a])) ++a;
        size_t b = a;
        auto isnum=[](char c){ return std::isdigit((unsigned char)c)||c=='-'||c=='+'||c=='.'||c=='e'||c=='E'; };
        while (b < src.size() && isnum(src[b])) ++b;
        try { return std::stod(src.substr(a, b-a)); } catch (...) { throw std::runtime_error("nan"); }
    };

    size_t pos = 0;
    int blocks = 0;
    while (true) {
        // Find next top-level DVCS_* key
        size_t kS = s.find('"', pos); if (kS == std::string::npos) break;
        size_t kE = s.find('"', kS+1); if (kE == std::string::npos) break;
        std::string key = s.substr(kS+1, kE-kS-1);
        pos = kE + 1;

        if (key.rfind("DVCS_", 0) != 0) continue;

        // Locate the `"data"` object for this block
        size_t dpos = s.find("\"data\"", kE);
        if (dpos == std::string::npos) fatal("Malformed cuts JSON near key: " + key);
        size_t dataObjStart = s.find('{', dpos);
        if (dataObjStart == std::string::npos) fatal("Malformed cuts JSON near key: " + key);

        // Extract the full braces-balanced `"data"` object
        int depth = 0;
        size_t i = dataObjStart;
        for (; i < s.size(); ++i) {
            if (s[i] == '{') depth++;
            else if (s[i] == '}') { depth--; if (!depth) { ++i; break; } }
        }
        if (depth != 0) fatal("Unbalanced braces in cuts JSON (data block for key " + key + ")");
        std::string data = s.substr(dataObjStart, i - dataObjStart);

        VarCutMap cuts;
        // Walk the top-level members of `"data"`: "<var>": { "mean": ..., "std": ... }
        size_t p = 0;
        while (true) {
            // find next quoted name
            size_t vS = data.find('"', p);
            if (vS == std::string::npos) break;
            size_t vE = data.find('"', vS+1);
            if (vE == std::string::npos) fatal("Malformed variable name in cuts JSON (key " + key + ")");
            std::string var = data.substr(vS+1, vE - vS - 1);

            // after name, expect ':' then an object { ... } for that variable
            size_t colon = data.find(':', vE+1);
            if (colon == std::string::npos) fatal("Malformed entry for var " + var + " (key " + key + ")");
            size_t objStart = data.find('{', colon+1);
            if (objStart == std::string::npos) {
                // If the next token isn’t an object, we likely hit an inner quoted token like "mean"
                // Advance and continue searching for the next top-level entry
                p = vE + 1;
                continue;
            }

            // Extract braces-balanced inner object for this variable
            int d2 = 0;
            size_t j = objStart;
            for (; j < data.size(); ++j) {
                if (data[j] == '{') d2++;
                else if (data[j] == '}') { d2--; if (!d2) { ++j; break; } }
            }
            if (d2 != 0) fatal("Unbalanced braces in var object for " + var + " (key " + key + ")");
            std::string varObj = data.substr(objStart, j - objStart);

            // Pull "mean" and "std" inside varObj
            size_t pm = varObj.find("\"mean\"");
            if (pm == std::string::npos) fatal("Missing 'mean' for var " + var + " (key " + key + ")");
            size_t cm = varObj.find(':', pm);
            if (cm == std::string::npos) fatal("Malformed 'mean' for var " + var + " (key " + key + ")");

            size_t ps = varObj.find("\"std\"");
            if (ps == std::string::npos) fatal("Missing 'std' for var " + var + " (key " + key + ")");
            size_t cs = varObj.find(':', ps);
            if (cs == std::string::npos) fatal("Malformed 'std' for var " + var + " (key " + key + ")");

            double mean = 0.0, stdev = 0.0;
            try { mean = parse_num_after_colon(varObj, cm); } catch(...) { fatal("Non-numeric 'mean' for var " + var + " (key " + key + ")"); }
            try { stdev = parse_num_after_colon(varObj, cs); } catch(...) { fatal("Non-numeric 'std' for var " + var + " (key " + key + ")"); }

            cuts[var] = Stats{mean, stdev};

            // advance after this variable's object
            p = j;
        }

        if (cuts.empty()) fatal("No variables parsed in cuts block for key " + key);
        out[key] = std::move(cuts);
        ++blocks;
    }

    if (blocks == 0) fatal("No DVCS_* blocks found in cuts JSON");
}

static inline bool within3Sigma(double v, const Stats& s){
    return (v >= s.mean - 3.0*s.std) && (v <= s.mean + 3.0*s.std);
}

static bool passes3SigmaCuts_STRICT(const VarCutMap& cuts,
                                    const std::map<std::string,double>& values,
                                    const std::string& period,
                                    const std::string& topoKey)
{
    for (const auto& kv : cuts) {
        auto it = values.find(kv.first);
        if (it == values.end()) {
            fatal("Cuts reference variable '" + kv.first + "' which is not available as a branch "
                  "(period=" + period + ", topo=" + topoKey + ")");
        }
        if (!within3Sigma(it->second, kv.second)) return false;
    }
    return true;
}

// ------------------------------------------------------------
// Branch binders (STRICT: require phi2, and branches used by cuts)
// ------------------------------------------------------------
struct BranchBinderDVCS {
    int detector1=0, detector2=0, helicity=0;
    double t1=0, open_angle_ep2=0, pTmiss=0;
    double x=0, Q2=0, phi2=0;

    // exclusivity-like quantities for 3σ on DVCS hypothesis
    // Provide all potential variables; if cuts reference any missing one -> FATAL in validation step
    double Emiss2=0, Mx2=0, Mx2_1=0, Mx2_2=0, theta_gamma_gamma=0, xF=0;

    std::set<std::string> boundD, boundI;

    void mustBindI(TTree* t, const char* n, int* a){ if(!t->GetBranch(n)) fatal(std::string("DVCS missing int branch '")+n+"'"); t->SetBranchAddress(n,a); boundI.insert(n); }
    void mustBindD(TTree* t, const char* n, double* a){ if(!t->GetBranch(n)) fatal(std::string("DVCS missing double branch '")+n+"'"); t->SetBranchAddress(n,a); boundD.insert(n); }

    void bind(TTree* t){
        if(!t) fatal("DVCS tree pointer is null");
        mustBindI(t,"detector1",&detector1);
        mustBindI(t,"detector2",&detector2);
        mustBindI(t,"helicity",&helicity);
        mustBindD(t,"t1",&t1);
        mustBindD(t,"open_angle_ep2",&open_angle_ep2);
        mustBindD(t,"pTmiss",&pTmiss);
        mustBindD(t,"x",&x);
        mustBindD(t,"Q2",&Q2);
        mustBindD(t,"phi2",&phi2);

        // Potential cut vars (must exist if referenced by cuts JSON)
        if (t->GetBranch("Emiss2"))            t->SetBranchAddress("Emiss2",&Emiss2), boundD.insert("Emiss2");
        if (t->GetBranch("Mx2"))               t->SetBranchAddress("Mx2",&Mx2), boundD.insert("Mx2");
        if (t->GetBranch("Mx2_1"))             t->SetBranchAddress("Mx2_1",&Mx2_1), boundD.insert("Mx2_1");
        if (t->GetBranch("Mx2_2"))             t->SetBranchAddress("Mx2_2",&Mx2_2), boundD.insert("Mx2_2");
        if (t->GetBranch("theta_gamma_gamma")) t->SetBranchAddress("theta_gamma_gamma",&theta_gamma_gamma), boundD.insert("theta_gamma_gamma");
        if (t->GetBranch("xF"))                t->SetBranchAddress("xF",&xF), boundD.insert("xF");
    }
    std::map<std::string,double> cutVals() const {
        std::map<std::string,double> v;
        // Always provide the core ones (φ cut is not applied by 3σ here)
        v["pTmiss"] = pTmiss;
        // Optional ones only if bound (validation will check if referenced)
        if (boundD.count("Emiss2")) v["Emiss2"]=Emiss2;
        if (boundD.count("Mx2")) v["Mx2"]=Mx2;
        if (boundD.count("Mx2_1")) v["Mx2_1"]=Mx2_1;
        if (boundD.count("Mx2_2")) v["Mx2_2"]=Mx2_2;
        if (boundD.count("theta_gamma_gamma")) v["theta_gamma_gamma"]=theta_gamma_gamma;
        if (boundD.count("xF")) v["xF"]=xF;
        return v;
    }
};

struct BranchBinderEPPI0Data { // has helicity
    int detector1=0, detector2=0, helicity=0;
    double t1=0, open_angle_ep2=0, pTmiss=0;
    double x=0, Q2=0, phi2=0;

    void mustBindI(TTree* t, const char* n, int* a){ if(!t->GetBranch(n)) fatal(std::string("eπ0 DATA missing int branch '")+n+"'"); t->SetBranchAddress(n,a); }
    void mustBindD(TTree* t, const char* n, double* a){ if(!t->GetBranch(n)) fatal(std::string("eπ0 DATA missing double branch '")+n+"'"); t->SetBranchAddress(n,a); }
    void bind(TTree* t){
        if(!t) fatal("eπ0 DATA tree pointer is null");
        mustBindI(t,"detector1",&detector1);
        mustBindI(t,"detector2",&detector2);
        mustBindI(t,"helicity",&helicity);
        mustBindD(t,"t1",&t1);
        mustBindD(t,"open_angle_ep2",&open_angle_ep2);
        mustBindD(t,"pTmiss",&pTmiss);
        mustBindD(t,"x",&x);
        mustBindD(t,"Q2",&Q2);
        mustBindD(t,"phi2",&phi2);
    }
};

struct BranchBinderEPPI0MC { // no helicity
    int detector1=0, detector2=0;
    double t1=0, open_angle_ep2=0, pTmiss=0;
    double x=0, Q2=0, phi2=0;

    // DVCS-hypothesis cuts variables that might be referenced
    double Emiss2=0, Mx2=0, Mx2_1=0, Mx2_2=0, theta_gamma_gamma=0, xF=0;
    std::set<std::string> boundD, boundI;

    void mustBindI(TTree* t, const char* n, int* a){ if(!t->GetBranch(n)) fatal(std::string("eπ0 MC missing int branch '")+n+"'"); t->SetBranchAddress(n,a); boundI.insert(n); }
    void mustBindD(TTree* t, const char* n, double* a){ if(!t->GetBranch(n)) fatal(std::string("eπ0 MC missing double branch '")+n+"'"); t->SetBranchAddress(n,a); boundD.insert(n); }
    void bind(TTree* t){
        if(!t) fatal("eπ0 MC tree pointer is null");
        mustBindI(t,"detector1",&detector1);
        mustBindI(t,"detector2",&detector2);
        mustBindD(t,"t1",&t1);
        mustBindD(t,"open_angle_ep2",&open_angle_ep2);
        mustBindD(t,"pTmiss",&pTmiss);
        mustBindD(t,"x",&x);
        mustBindD(t,"Q2",&Q2);
        mustBindD(t,"phi2",&phi2);

        // Potential DVCS-hypothesis cut vars
        if (t->GetBranch("Emiss2"))            t->SetBranchAddress("Emiss2",&Emiss2), boundD.insert("Emiss2");
        if (t->GetBranch("Mx2"))               t->SetBranchAddress("Mx2",&Mx2), boundD.insert("Mx2");
        if (t->GetBranch("Mx2_1"))             t->SetBranchAddress("Mx2_1",&Mx2_1), boundD.insert("Mx2_1");
        if (t->GetBranch("Mx2_2"))             t->SetBranchAddress("Mx2_2",&Mx2_2), boundD.insert("Mx2_2");
        if (t->GetBranch("theta_gamma_gamma")) t->SetBranchAddress("theta_gamma_gamma",&theta_gamma_gamma), boundD.insert("theta_gamma_gamma");
        if (t->GetBranch("xF"))                t->SetBranchAddress("xF",&xF), boundD.insert("xF");
    }
    std::map<std::string,double> cutValsForDVCS() const {
        std::map<std::string,double> v;
        v["pTmiss"] = pTmiss;
        if (boundD.count("Emiss2")) v["Emiss2"]=Emiss2;
        if (boundD.count("Mx2")) v["Mx2"]=Mx2;
        if (boundD.count("Mx2_1")) v["Mx2_1"]=Mx2_1;
        if (boundD.count("Mx2_2")) v["Mx2_2"]=Mx2_2;
        if (boundD.count("theta_gamma_gamma")) v["theta_gamma_gamma"]=theta_gamma_gamma;
        if (boundD.count("xF")) v["xF"]=xF;
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
// Writers (STRICT)
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
    if (!ofs) fatal(std::string("open for write failed: ") + out_path);
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
    ofs.close();
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
    if (!ofs) fatal(std::string("open for write failed: ") + out_path);
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
    ofs.close();
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

// ------------------------------------------------------------
// Plotting
// ------------------------------------------------------------
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
    namespace fs = std::filesystem;

    std::error_code ec;
    fs::create_directories(out_dir, ec);
    if (ec) {
        fatal(std::string("[pi0_contam][FATAL] Cannot create directory: ") + out_dir +
              " (" + ec.message() + ")");
    }

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
        head.DrawLatex(0.5, 0.55, Form("%s   x_{B} #in [%.3g, %.3g]", name.c_str(),
                                       xB_bins[ix].first, xB_bins[ix].second));

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

        // Save and verify on disk (fail-fast)
        const std::string fpath = out_dir + "/plot_contamination_" + name + "_xB_" + std::to_string(ix) + ".png";
        c->SaveAs(fpath.c_str());

        std::error_code fec;
        const bool exists = fs::exists(fpath, fec);
        const auto size   = exists ? fs::file_size(fpath, fec) : 0ULL;
        if (!exists || size == 0 || fec) {
            delete c;
            std::ostringstream em;
            em << "[pi0_contam][FATAL] Failed to save plot: " << fpath
               << " (exists=" << std::boolalpha << exists
               << ", size=" << size
               << ", ec=" << (fec ? fec.message() : "ok") << ")";
            fatal(em.str());
        }

        delete c;
    }
}

} // namespace

// ------------------------------------------------------------
// Core (STRICT)
// ------------------------------------------------------------
void compute_pi0_contamination_helicity(
    const std::vector<std::string>& periods,
    const std::vector<std::string>& topologies,
    const std::vector<Binning>& binning_scheme,
    const std::map<std::string, TTree*>& dvcsDataTrees,   // key: DVCS_<...>
    const std::map<std::string, TTree*>& eppi0DataTrees,  // key: DVCS_<...>_eppi0
    const std::map<std::string, TTree*>& eppi0RecMcTrees, // key: DVCS_<...>_rec_mc
    const std::map<std::string, TTree*>& eppi0BkgTrees,   // key: DVCS_<...>_bkg
    const std::string& combined_cuts_json,
    const std::string& out_root_dir)
{
    namespace fs = std::filesystem;

    if (topologies.empty()) fatal("Topologies list is empty");
    const auto xB_bins = uniqueRanges(binning_scheme,'x');
    const auto Q2_bins = uniqueRanges(binning_scheme,'Q');
    const auto t_bins  = uniqueRanges(binning_scheme,'t');
    if (xB_bins.empty() || Q2_bins.empty() || t_bins.empty())
        fatal("Empty bin ranges (xB/Q2/t). Check binning_scheme.");

    // Load 3σ cuts (STRICT)
    PeriodTopoCuts cuts;
    loadCombinedCuts_STRICT(combined_cuts_json, cuts);

    // EXACT key: if runTag already starts with "DVCS_", do not add; no titlecasing.
    auto dvcsCutsKey = [&](const std::string& runTag, const std::string& topoKey)->std::string {
        std::string base = runTag;
        if (base.rfind("DVCS_", 0) != 0) base = std::string("DVCS_") + base;
        return base + "_" + topoKey;
    };

    // Output paths (STRICT)
    const fs::path root(out_root_dir);
    const fs::path jsons_dir    = root / "jsons" / "contamination";
    const fs::path combined_out = root / "jsons" / "pi0_contamination_combined.json";
    const fs::path plots_root   = root / "contamination_plots";
    std::error_code ec;
    if (!fs::create_directories(jsons_dir, ec) && ec)
        fatal(std::string("Cannot create jsons/contamination dir: ")+jsons_dir.string()+" ("+ec.message()+")");
    ec.clear();
    if (!fs::create_directories(plots_root, ec) && ec)
        fatal(std::string("Cannot create plots root: ")+plots_root.string()+" ("+ec.message()+")");

    // Keep all per-period tables for combined build and plots
    std::map<std::string, std::map<BinKey, BinCounts>> groupTables;

    // ---------- Per-period loops (STRICT) ----------
    for (const auto& period : periods) {
        const std::string runTag = periodToRunTagKey(period); // should already be canonical (e.g. DVCS_Sp18_inb)
        const std::string key_dvcs     = runTag;               // DVCS data
        const std::string key_pi0_data = runTag + "_eppi0";    // eπ0 data
        const std::string key_pi0_reco = runTag + "_rec_mc";   // eπ0 reco MC
        const std::string key_pi0_bkg  = runTag + "_bkg";      // eπ0→DVCS bkg MC

        auto itD  = dvcsDataTrees.find(key_dvcs);
        auto itED = eppi0DataTrees.find(key_pi0_data);
        auto itRM = eppi0RecMcTrees.find(key_pi0_reco);
        auto itBM = eppi0BkgTrees.find(key_pi0_bkg);

        if (itD==dvcsDataTrees.end() || !itD->second)
            fatal("Missing DVCS DATA for " + period + " (key: " + key_dvcs + ")");
        if (itED==eppi0DataTrees.end() || !itED->second)
            fatal("Missing eπ0 DATA for " + period + " (key: " + key_pi0_data + ")");
        if (itRM==eppi0RecMcTrees.end() || !itRM->second)
            fatal("Missing eπ0 RECO MC for " + period + " (key: " + key_pi0_reco + ")");
        if (itBM==eppi0BkgTrees.end() || !itBM->second)
            fatal("Missing eπ0→DVCS BKG MC for " + period + " (key: " + key_pi0_bkg + ")");

        TTree* t_dvcs     = itD->second;
        TTree* t_pi0_data = itED->second;
        TTree* t_pi0_reco = itRM->second;
        TTree* t_pi0_bkg  = itBM->second;

        std::map<BinKey, BinCounts> table;

        // DVCS data → N_data^±
        {
            BranchBinderDVCS b; b.bind(t_dvcs);

            // Validate that any cut variable referenced by JSON exists in DVCS tree
            for (const auto& topoStr : topologies) {
                const std::string topoKey = topoToKey(topoStr);
                const std::string cKey    = dvcsCutsKey(runTag, topoKey);
                auto itC = cuts.find(cKey);
                if (itC == cuts.end())
                    fatal("Cuts block not found for key: " + cKey + " (period: " + period + ")");
                for (const auto& v : itC->second) {
                    const std::string& var = v.first;
                    // Branch must be bound if referenced
                    if (b.boundD.count(var)==0) {
                        fatal("Cuts reference variable '" + var + "' but DVCS tree does not provide it "
                              "(period=" + period + ", topo=" + topoKey + ")");
                    }
                }
            }

            const Long64_t nent = t_dvcs->GetEntries();
            for (Long64_t i=0;i<nent;++i) {
                t_dvcs->GetEntry(i);
                if (!applyKinematicCuts_simple(b.t1,b.open_angle_ep2,b.pTmiss)) continue;
                if (b.helicity!=+1 && b.helicity!=-1) continue;

                std::string usedTopoKey;
                bool matched=false;
                for (const auto& topoStr : topologies) {
                    if (passesTopology_simple(b.detector1,b.detector2,topoStr)) {
                        usedTopoKey = topoToKey(topoStr); matched=true; break;
                    }
                }
                if (!matched) continue;

                const std::string cKey = dvcsCutsKey(runTag, usedTopoKey);
                auto itC = cuts.find(cKey);
                if (itC == cuts.end()) fatal("Missing cuts block (post-validation) for key " + cKey);
                if (!passes3SigmaCuts_STRICT(itC->second, b.cutVals(), period, usedTopoKey)) continue;

                double xB=b.x, Q2=b.Q2, tt=std::fabs(b.t1), phi=b.phi2;
                if (!std::isfinite(xB)||!std::isfinite(Q2)||!std::isfinite(tt)||!std::isfinite(phi)) continue;
                int ix=findBin(xB,xB_bins), iQ=findBin(Q2,Q2_bins), itb=findBin(tt,t_bins), ip=phiToBin(phi);
                if (ix<0||iQ<0||itb<0||ip<0||ip>=N_PHI_BINS) continue;

                BinKey key(ix,iQ,itb,ip);
                if (b.helicity==+1) table[key].N_data.plus++; else table[key].N_data.minus++;
            }
        }

        // eπ0→DVCS mis-ID MC → N_pi0_mc
        {
            BranchBinderEPPI0MC b; b.bind(t_pi0_bkg);

            // Validate referenced DVCS cuts variables also exist on this MC if cuts applied to MC
            for (const auto& topoStr : topologies) {
                const std::string topoKey = topoToKey(topoStr);
                const std::string cKey    = dvcsCutsKey(runTag, topoKey);
                auto itC = cuts.find(cKey);
                if (itC == cuts.end()) fatal("Cuts block not found for key: " + cKey);
                for (const auto& v : itC->second) {
                    const std::string& var = v.first;
                    if (var=="pTmiss") continue; // guaranteed
                    if (b.boundD.count(var)==0) {
                        fatal("Cuts reference variable '" + var + "' but eπ0_bkg MC does not provide it "
                              "(period=" + period + ", topo=" + topoKey + ")");
                    }
                }
            }

            const Long64_t nent = t_pi0_bkg->GetEntries();
            for (Long64_t i=0;i<nent;++i) {
                t_pi0_bkg->GetEntry(i);
                if (!applyKinematicCuts_simple(b.t1,b.open_angle_ep2,b.pTmiss)) continue;

                bool match=false; std::string usedTopoKey;
                for (const auto& topoStr : topologies) {
                    if (passesTopology_simple(b.detector1,b.detector2,topoStr)) { match=true; usedTopoKey = topoToKey(topoStr); break; }
                }
                if (!match) continue;

                const std::string cKey = dvcsCutsKey(runTag, usedTopoKey);
                auto itC = cuts.find(cKey);
                if (itC == cuts.end()) fatal("Missing cuts block (post-validation) for key " + cKey);
                if (!passes3SigmaCuts_STRICT(itC->second, b.cutValsForDVCS(), period, usedTopoKey)) continue;

                double xB=b.x, Q2=b.Q2, tt=std::fabs(b.t1), phi=b.phi2;
                if (!std::isfinite(xB)||!std::isfinite(Q2)||!std::isfinite(tt)||!std::isfinite(phi)) continue;
                int ix=findBin(xB,xB_bins), iQ=findBin(Q2,Q2_bins), itb=findBin(tt,t_bins), ip=phiToBin(phi);
                if (ix<0||iQ<0||itb<0||ip<0||ip>=N_PHI_BINS) continue;

                table[BinKey(ix,iQ,itb,ip)].N_pi0_mc++;
            }
        }

        // eπ0 reco MC → N_pi0_reco
        {
            BranchBinderEPPI0MC b; b.bind(t_pi0_reco);
            const Long64_t nent = t_pi0_reco->GetEntries();
            for (Long64_t i=0;i<nent;++i) {
                t_pi0_reco->GetEntry(i);
                if (!applyKinematicCuts_simple(b.t1,b.open_angle_ep2,b.pTmiss)) continue;

                bool match=false;
                for (const auto& topoStr : topologies) {
                    if (passesTopology_simple(b.detector1,b.detector2,topoStr)) { match=true; break; }
                }
                if (!match) continue;

                double xB=b.x, Q2=b.Q2, tt=std::fabs(b.t1), phi=b.phi2;
                if (!std::isfinite(xB)||!std::isfinite(Q2)||!std::isfinite(tt)||!std::isfinite(phi)) continue;
                int ix=findBin(xB,xB_bins), iQ=findBin(Q2,Q2_bins), itb=findBin(tt,t_bins), ip=phiToBin(phi);
                if (ix<0||iQ<0||itb<0||ip<0||ip>=N_PHI_BINS) continue;

                table[BinKey(ix,iQ,itb,ip)].N_pi0_reco++;
            }
        }

        // eπ0 DATA → N_pi0_exp^±
        {
            BranchBinderEPPI0Data b; b.bind(t_pi0_data);
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

                double xB=b.x, Q2=b.Q2, tt=std::fabs(b.t1), phi=b.phi2;
                if (!std::isfinite(xB)||!std::isfinite(Q2)||!std::isfinite(tt)||!std::isfinite(phi)) continue;
                int ix=findBin(xB,xB_bins), iQ=findBin(Q2,Q2_bins), itb=findBin(tt,t_bins), ip=phiToBin(phi);
                if (ix<0||iQ<0||itb<0||ip<0||ip>=N_PHI_BINS) continue;

                BinKey key(ix,iQ,itb,ip);
                if (b.helicity==+1) table[key].N_pi0_exp.plus++; else table[key].N_pi0_exp.minus++;
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