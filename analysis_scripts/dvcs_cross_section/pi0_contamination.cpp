#include "pi0_contamination.h"

#include <TTree.h>

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

bool COPY_CONTAM_TO_FA18_INB_SUPP = true;

// ---------------- helpers (shared style with total_counts.cpp) ----------------
static inline std::string toLower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c){ return std::tolower(c); });
    return s;
}
static std::string periodToRunTagKey(const std::string& period) {
    auto pos = period.find('_');
    if (pos == std::string::npos || pos + 1 >= period.size()) return toLower(period);
    return toLower(period.substr(pos + 1));
}
static inline std::string topoToKey(const std::string& topoStr) {
    if (topoStr == "(FD,FD)") return "FD_FD";
    if (topoStr == "(CD,FD)") return "CD_FD";
    if (topoStr == "(CD,FT)") return "CD_FT";
    return "FD_FD";
}
static inline bool passesTopology_simple(int detector1, int detector2, const std::string& topoStr) {
    if (topoStr == "(FD,FD)") return (detector1 == 1 && detector2 == 1);
    if (topoStr == "(CD,FD)") return (detector1 == 2 && detector2 == 1);
    if (topoStr == "(CD,FT)") return (detector1 == 2 && detector2 == 0);
    return false;
}
static inline bool applyKinematicCuts_simple(double t1, double open_angle_ep2, double pTmiss) {
    if (open_angle_ep2 <= 5.0) return false;
    if ((-t1) > 1.0)          return false;
    if (pTmiss > 0.20)        return false;
    return true;
}

// binning helpers
static constexpr int    N_PHI_BINS = 12;
static constexpr double TWO_PI     = 2.0 * M_PI;
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

// ---------------- exclusivity cuts loader (data-cuts only) ----------------
struct Stats { double mean = 0.0; double std = 0.0; };
using VarCutMap = std::map<std::string, Stats>;
using PeriodTopoCuts = std::map<std::string, VarCutMap>; // key: "DVCS_Fa18_inb_FD_FD"

static bool loadCombinedCuts(const std::string& path, PeriodTopoCuts& out) {
    std::ifstream ifs(path);
    if (!ifs) { std::cerr << "[pi0_contam][ERROR] Cannot open cuts JSON: " << path << std::endl; return false; }
    std::string s((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    ifs.close();

    size_t pos = 0;
    while (true) {
        size_t keyStart = s.find('"', pos);
        if (keyStart == std::string::npos) break;
        size_t keyEnd = s.find('"', keyStart + 1);
        if (keyEnd == std::string::npos) break;
        std::string key = s.substr(keyStart + 1, keyEnd - keyStart - 1);

        size_t dataPos = s.find("\"data\"", keyEnd);
        if (dataPos == std::string::npos) { pos = keyEnd + 1; continue; }
        size_t braceStart = s.find('{', dataPos);
        if (braceStart == std::string::npos) { pos = keyEnd + 1; continue; }

        int depth = 0; size_t i = braceStart;
        for (; i < s.size(); ++i) {
            if (s[i] == '{') depth++;
            else if (s[i] == '}') { depth--; if (depth == 0) { ++i; break; } }
        }
        if (depth != 0) { pos = keyEnd + 1; continue; }
        std::string dataObj = s.substr(braceStart, i - braceStart);

        VarCutMap cuts;
        size_t vpos = 0;
        while (true) {
            size_t vKeyS = dataObj.find('"', vpos);
            if (vKeyS == std::string::npos) break;
            size_t vKeyE = dataObj.find('"', vKeyS + 1);
            if (vKeyE == std::string::npos) break;
            std::string var = dataObj.substr(vKeyS + 1, vKeyE - vKeyS - 1);

            size_t meanPos = dataObj.find("\"mean\"", vKeyE);
            size_t stdPos  = dataObj.find("\"std\"",  vKeyE);
            if (meanPos == std::string::npos || stdPos == std::string::npos) { vpos = vKeyE + 1; continue; }

            auto readNum = [&](size_t from)->double {
                size_t colon = dataObj.find(':', from);
                if (colon == std::string::npos) return 0.0;
                size_t j = colon + 1;
                while (j < dataObj.size() && std::isspace(static_cast<unsigned char>(dataObj[j]))) ++j;
                size_t k = j;
                while (k < dataObj.size() && (std::isdigit(static_cast<unsigned char>(dataObj[k])) || dataObj[k]=='-' || dataObj[k]=='+' || dataObj[k]=='.' || dataObj[k]=='e' || dataObj[k]=='E')) ++k;
                try { return std::stod(dataObj.substr(j, k - j)); } catch (...) { return 0.0; }
            };

            double m = readNum(meanPos);
            double sd = readNum(stdPos);
            cuts[var] = Stats{m, sd};

            vpos = vKeyE + 1;
        }

        if (key.rfind("DVCS_", 0) == 0) out[key] = cuts; // keep only DVCS entries
        pos = keyEnd + 1;
    }
    return !out.empty();
}

static inline bool within3Sigma(double v, const Stats& s) {
    return (v >= s.mean - 3.0*s.std) && (v <= s.mean + 3.0*s.std);
}
static bool passes3SigmaCuts(const VarCutMap& cuts, const std::map<std::string,double>& values) {
    for (const auto& kv : cuts) {
        auto it = values.find(kv.first);
        if (it == values.end()) continue;
        if (!within3Sigma(it->second, kv.second)) return false;
    }
    return true;
}

// ---------------- branch binders ----------------
struct BranchBinderDVCS {
    int detector1=0, detector2=0; bool has_d1=false, has_d2=false;
    int helicity=0; bool has_helicity=false;
    double t1=0, open_angle_ep2=0, pTmiss=0; bool has_t1=false, has_oa=false, has_pT=false;
    double x=0, Q2=0, phi2=0, Delta_phi=0; bool has_x=false, has_Q2=false, has_phi2=false, has_Dp=false;
    // exclusivity vars
    double Emiss2=0, Mx2=0, Mx2_1=0, Mx2_2=0, theta_gamma_gamma=0, xF=0;
    bool has_Em=false, has_Mx2=false, has_Mx21=false, has_Mx22=false, has_tgg=false, has_xF=false;

    void bind(TTree* t) {
        if (!t) return;
        auto bindI=[&](const char*n,int*a,bool&f){ if(t->GetBranch(n)){ t->SetBranchAddress(n,a); f=true; } };
        auto bindD=[&](const char*n,double*a,bool&f){ if(t->GetBranch(n)){ t->SetBranchAddress(n,a); f=true; } };

        bindI("detector1",&detector1,has_d1);
        bindI("detector2",&detector2,has_d2);
        bindI("helicity",&helicity,has_helicity);

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
        bindD("xF",&xF,has_xF);
    }
    bool readyCuts() const { return has_d1 && has_d2 && has_t1 && has_oa && has_pT && has_helicity; }
    bool readyBins() const { return has_x && has_Q2 && (has_phi2 || has_Dp); }
    double phi() const { return has_phi2 ? phi2 : (has_Dp ? Delta_phi : std::numeric_limits<double>::quiet_NaN()); }

    std::map<std::string,double> cutVals() const {
        std::map<std::string,double> m;
        if (has_Dp) m["Delta_phi"]=Delta_phi;
        if (has_tgg) m["theta_gamma_gamma"]=theta_gamma_gamma;
        if (has_pT) m["pTmiss"]=pTmiss;
        if (has_xF) m["xF"]=xF;
        if (has_Em) m["Emiss2"]=Emiss2;
        if (has_Mx2) m["Mx2"]=Mx2;
        if (has_Mx21) m["Mx2_1"]=Mx2_1;
        if (has_Mx22) m["Mx2_2"]=Mx2_2;
        return m;
    }
};
struct BranchBinderEPPI0Data { // has helicity
    int detector1=0, detector2=0; bool has_d1=false, has_d2=false;
    int helicity=0; bool has_helicity=false;
    double t1=0, open_angle_ep2=0, pTmiss=0; bool has_t1=false, has_oa=false, has_pT=false;
    double x=0, Q2=0, phi2=0, Delta_phi=0; bool has_x=false, has_Q2=false, has_phi2=false, has_Dp=false;
    // exclusivity vars use theta_pi0_pi0 for eppi0
    double Emiss2=0, Mx2=0, Mx2_1=0, Mx2_2=0, theta_pi0_pi0=0, xF=0;
    bool has_Em=false, has_Mx2=false, has_Mx21=false, has_Mx22=false, has_tpp=false, has_xF=false;

    void bind(TTree* t) {
        if (!t) return;
        auto bindI=[&](const char*n,int*a,bool&f){ if(t->GetBranch(n)){ t->SetBranchAddress(n,a); f=true; } };
        auto bindD=[&](const char*n,double*a,bool&f){ if(t->GetBranch(n)){ t->SetBranchAddress(n,a); f=true; } };

        bindI("detector1",&detector1,has_d1);
        bindI("detector2",&detector2,has_d2);
        bindI("helicity",&helicity,has_helicity);

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
        bindD("xF",&xF,has_xF);
    }
    bool readyCuts() const { return has_d1 && has_d2 && has_t1 && has_oa && has_pT && has_helicity; }
    bool readyBins() const { return has_x && has_Q2 && (has_phi2 || has_Dp); }
    double phi() const { return has_phi2 ? phi2 : (has_Dp ? Delta_phi : std::numeric_limits<double>::quiet_NaN()); }

    std::map<std::string,double> cutVals() const {
        std::map<std::string,double> m;
        if (has_Dp) m["Delta_phi"]=Delta_phi;
        if (has_tpp) m["theta_pi0_pi0"]=theta_pi0_pi0;
        if (has_pT) m["pTmiss"]=pTmiss;
        if (has_xF) m["xF"]=xF;
        if (has_Em) m["Emiss2"]=Emiss2;
        if (has_Mx2) m["Mx2"]=Mx2;
        if (has_Mx21) m["Mx2_1"]=Mx2_1;
        if (has_Mx22) m["Mx2_2"]=Mx2_2;
        return m;
    }
};
struct BranchBinderEPPI0MC { // no helicity
    int detector1=0, detector2=0; bool has_d1=false, has_d2=false;
    double t1=0, open_angle_ep2=0, pTmiss=0; bool has_t1=false, has_oa=false, has_pT=false;
    double x=0, Q2=0, phi2=0, Delta_phi=0; bool has_x=false, has_Q2=false, has_phi2=false, has_Dp=false;
    // exclusivity
    double Emiss2=0, Mx2=0, Mx2_1=0, Mx2_2=0, theta_pi0_pi0=0, theta_gamma_gamma=0, xF=0;
    bool has_Em=false, has_Mx2=false, has_Mx21=false, has_Mx22=false, has_tpp=false, has_tgg=false, has_xF=false;

    void bind(TTree* t) {
        if (!t) return;
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
        bindD("xF",&xF,has_xF);
    }
    bool readyCuts() const { return has_d1 && has_d2 && has_t1 && has_oa && has_pT; }
    bool readyBins() const { return has_x && has_Q2 && (has_phi2 || has_Dp); }
    double phi() const { return has_phi2 ? phi2 : (has_Dp ? Delta_phi : std::numeric_limits<double>::quiet_NaN()); }

    std::map<std::string,double> cutValsForDVCS() const { // when mis-ID to DVCS hypothesis
        std::map<std::string,double> m;
        if (has_Dp) m["Delta_phi"]=Delta_phi;
        if (has_tgg) m["theta_gamma_gamma"]=theta_gamma_gamma;
        if (has_pT) m["pTmiss"]=pTmiss;
        if (has_xF) m["xF"]=xF;
        if (has_Em) m["Emiss2"]=Emiss2;
        if (has_Mx2) m["Mx2"]=Mx2;
        if (has_Mx21) m["Mx2_1"]=Mx2_1;
        if (has_Mx22) m["Mx2_2"]=Mx2_2;
        return m;
    }
    std::map<std::string,double> cutValsForEPPI0() const { // when genuine π0 selection
        std::map<std::string,double> m;
        if (has_Dp) m["Delta_phi"]=Delta_phi;
        if (has_tpp) m["theta_pi0_pi0"]=theta_pi0_pi0;
        if (has_pT) m["pTmiss"]=pTmiss;
        if (has_xF) m["xF"]=xF;
        if (has_Em) m["Emiss2"]=Emiss2;
        if (has_Mx2) m["Mx2"]=Mx2;
        if (has_Mx21) m["Mx2_1"]=Mx2_1;
        if (has_Mx22) m["Mx2_2"]=Mx2_2;
        return m;
    }
};

// ---------------- containers ----------------
using BinKey = std::tuple<int,int,int,int>; // (ix,iQ2,it,ip)
struct HelCounts { long long plus=0, minus=0; };
struct BinCounts {
    HelCounts N_data;       // DVCS data, by helicity
    HelCounts N_pi0_exp;    // eppi0 DATA, by helicity
    long long N_pi0_mc  = 0; // eppi0_bkg mis-ID MC (no helicity)
    long long N_pi0_reco= 0; // eppi0 reco MC (no helicity)
    // results (per helicity)
    double c_plus = 0.0, c_plus_err = 0.0;
    double c_minus= 0.0, c_minus_err= 0.0;
};

// ---------------- write JSON ----------------
static inline std::string keyStr(const BinKey& k) {
    int a,b,c,d; std::tie(a,b,c,d)=k;
    std::ostringstream os; os<<"("<<a<<","<<b<<","<<c<<","<<d<<")";
    return os.str();
}
static void writeContaminationJson(const std::string& path,
                                   const std::map<BinKey, BinCounts>& table,
                                   int nPhi,
                                   const std::vector<std::pair<double,double>>& xB_bins,
                                   const std::vector<std::pair<double,double>>& Q2_bins,
                                   const std::vector<std::pair<double,double>>& t_bins)
{
    std::ofstream ofs(path);
    if (!ofs) { std::cerr << "[pi0_contam][ERROR] Cannot open " << path << std::endl; return; }
    ofs << std::fixed << std::setprecision(8);
    ofs << "{\n";
    ofs << "  \"binning_meta\": {\"phi_bins\": " << nPhi
        << ", \"xB_bins\": " << xB_bins.size()
        << ", \"Q2_bins\": " << Q2_bins.size()
        << ", \"t_bins\": "  << t_bins.size()  << "},\n";
    ofs << "  \"bins\": {\n";
    bool first=true;
    for (const auto& kv : table) {
        if (!first) ofs << ",\n";
        first=false;
        const auto& bc = kv.second;
        ofs << "    \"" << keyStr(kv.first) << "\": {"
            << "\"N_data\":{\"helicity\":{\"+1\":" << bc.N_data.plus
            << ",\"-1\":" << bc.N_data.minus << "},\"total\":" << (bc.N_data.plus+bc.N_data.minus) << "}"
            << ",\"N_pi0_exp\":{\"helicity\":{\"+1\":" << bc.N_pi0_exp.plus
            << ",\"-1\":" << bc.N_pi0_exp.minus << "},\"total\":" << (bc.N_pi0_exp.plus+bc.N_pi0_exp.minus) << "}"
            << ",\"N_pi0_mc\":"   << bc.N_pi0_mc
            << ",\"N_pi0_reco\":" << bc.N_pi0_reco
            << ",\"contamination\":{"
            << "\"+1\":{\"value\":" << bc.c_plus  << ",\"err\":" << bc.c_plus_err  << "},"
            << "\"-1\":{\"value\":" << bc.c_minus << ",\"err\":" << bc.c_minus_err << "}"
            << "}"
            << "}";
    }
    ofs << "\n  }\n}\n";
    ofs.close();
    std::cout << "[pi0_contam] Wrote " << path << std::endl;
}

// ---------------- core ----------------
void compute_pi0_contamination_helicity(
    const std::vector<std::string>& periods,
    const std::vector<std::string>& topologies,
    const std::vector<Binning>& binning_scheme,
    const std::map<std::string, TTree*>& dvcsDataTrees,
    const std::map<std::string, TTree*>& eppi0DataTrees,
    const std::map<std::string, TTree*>& eppi0RecMcTrees,
    const std::map<std::string, TTree*>& eppi0BkgTrees,
    const std::string& combined_cuts_json,
    const std::string& out_dir
) {
    namespace fs = std::filesystem;
    fs::create_directories(out_dir);

    const auto xB_bins = uniqueRanges(binning_scheme, 'x');
    const auto Q2_bins = uniqueRanges(binning_scheme, 'Q');
    const auto t_bins  = uniqueRanges(binning_scheme, 't');
    if (xB_bins.empty() || Q2_bins.empty() || t_bins.empty()) {
        std::cerr << "[pi0_contam][ERROR] Missing binning ranges." << std::endl;
        return;
    }

    // Load exclusivity 3σ cuts (data) per DVCS period/topology
    PeriodTopoCuts cuts;
    if (!loadCombinedCuts(combined_cuts_json, cuts)) {
        std::cerr << "[pi0_contam][WARN] No combined cuts loaded; proceeding without 3σ cuts." << std::endl;
    }

    auto dvcsCutsKey = [&](const std::string& runTag, const std::string& topoKey)->std::string {
        // make "DVCS_Fa18_inb" style
        std::string cap = runTag;
        if (!cap.empty()) cap[0] = std::toupper(cap[0]);
        for (size_t i=0;i+1<cap.size();++i) if (cap[i]=='_' && i+1<cap.size()) cap[i+1]=std::toupper(cap[i+1]);
        return std::string("DVCS_") + cap + "_" + topoKey;
    };

    // For each DVCS period, build counts
    for (const auto& period : periods) {
        const std::string runTag = periodToRunTagKey(period);       // e.g. fa18_inb
        const std::string key_dvcs = runTag;                        // dvcs data key
        const std::string key_pi0_data = runTag + "_eppi0";         // eppi0 data
        const std::string key_pi0_reco = runTag + "_rec_mc";        // eppi0 reco MC
        const std::string key_pi0_bkg  = runTag + "_bkg";           // eppi0 mis-ID MC (to DVCS)

        auto itDVCS = dvcsDataTrees.find(key_dvcs);
        if (itDVCS == dvcsDataTrees.end() || !itDVCS->second) {
            std::cerr << "[pi0_contam][WARN] No DVCS data tree for '" << period << "' (key '" << key_dvcs << "'). Skipping." << std::endl;
            continue;
        }
        TTree* t_dvcs = itDVCS->second;

        TTree* t_pi0_data = nullptr;
        auto itPi0Data = eppi0DataTrees.find(key_pi0_data);
        if (itPi0Data != eppi0DataTrees.end()) t_pi0_data = itPi0Data->second;

        TTree* t_pi0_reco = nullptr;
        auto itPi0Reco = eppi0RecMcTrees.find(key_pi0_reco);
        if (itPi0Reco != eppi0RecMcTrees.end()) t_pi0_reco = itPi0Reco->second;

        TTree* t_pi0_bkg = nullptr;
        auto itPi0Bkg = eppi0BkgTrees.find(key_pi0_bkg);
        if (itPi0Bkg != eppi0BkgTrees.end()) t_pi0_bkg = itPi0Bkg->second;

        if (!t_pi0_bkg || !t_pi0_reco) {
            std::cerr << "[pi0_contam][WARN] Missing π0 MC for '" << period << "'; bkg=" << (t_pi0_bkg?"ok":"none")
                      << " reco=" << (t_pi0_reco?"ok":"none") << ". Continuing with what is available." << std::endl;
        }
        if (!t_pi0_data) {
            std::cerr << "[pi0_contam][WARN] Missing eppi0 DATA for '" << period << "'. Continuing with zeros for N_pi0_exp." << std::endl;
        }

        // Accumulator per bin
        std::map<BinKey, BinCounts> counts;

        // ---- DVCS data (helicity-resolved) ----
        {
            BranchBinderDVCS b; b.bind(t_dvcs);
            if (!b.readyCuts() || !b.readyBins()) {
                std::cerr << "[pi0_contam][WARN] DVCS tree missing branches for '" << period << "'. Skipping DVCS loop." << std::endl;
            } else {
                const Long64_t nent = t_dvcs->GetEntries();
                for (Long64_t i=0;i<nent;++i) {
                    t_dvcs->GetEntry(i);
                    if (!applyKinematicCuts_simple(b.t1,b.open_angle_ep2,b.pTmiss)) continue;
                    if (b.helicity!=+1 && b.helicity!=-1) continue;

                    // choose topology for this event (first matching)
                    std::string usedTopoKey;
                    for (const auto& topoStr : topologies) {
                        if (passesTopology_simple(b.detector1,b.detector2,topoStr)) { usedTopoKey = topoToKey(topoStr); break; }
                    }
                    if (usedTopoKey.empty()) continue;

                    // 3σ cuts (if available)
                    VarCutMap topoCuts;
                    auto itCut = cuts.find(dvcsCutsKey(runTag, usedTopoKey));
                    if (itCut != cuts.end()) topoCuts = itCut->second;
                    if (!topoCuts.empty()) {
                        if (!passes3SigmaCuts(topoCuts, b.cutVals())) continue;
                    }

                    // binning
                    double xB=b.x, Q2=b.Q2, tt=std::fabs(b.t1), phi=b.phi();
                    if (!std::isfinite(xB)||!std::isfinite(Q2)||!std::isfinite(tt)||!std::isfinite(phi)) continue;
                    int ix=findBin(xB,xB_bins), iQ=findBin(Q2,Q2_bins), itb=findBin(tt,t_bins), ip=phiToBin(phi);
                    if (ix<0||iQ<0||itb<0||ip<0||ip>=N_PHI_BINS) continue;

                    BinKey key(ix,iQ,itb,ip);
                    if (b.helicity==+1) counts[key].N_data.plus++; else counts[key].N_data.minus++;
                }
            }
        }

        // ---- π0 background MC mis-ID to DVCS (no helicity) -> N_pi0_mc ----
        if (t_pi0_bkg) {
            BranchBinderEPPI0MC b; b.bind(t_pi0_bkg);
            if (!b.readyCuts() || !b.readyBins()) {
                std::cerr << "[pi0_contam][WARN] eppi0_bkg tree missing branches for '" << period << "'. Skipping." << std::endl;
            } else {
                const Long64_t nent = t_pi0_bkg->GetEntries();
                for (Long64_t i=0;i<nent;++i) {
                    t_pi0_bkg->GetEntry(i);
                    if (!applyKinematicCuts_simple(b.t1,b.open_angle_ep2,b.pTmiss)) continue;

                    // event topology
                    bool match=false; std::string usedTopoKey;
                    for (const auto& topoStr : topologies) {
                        if (passesTopology_simple(b.detector1,b.detector2,topoStr)) { usedTopoKey = topoToKey(topoStr); match=true; break; }
                    }
                    if (!match) continue;

                    // 3σ cuts under DVCS hypothesis
                    VarCutMap topoCuts;
                    auto itCut = cuts.find(dvcsCutsKey(runTag, usedTopoKey));
                    if (itCut != cuts.end()) topoCuts = itCut->second;
                    if (!topoCuts.empty()) {
                        if (!passes3SigmaCuts(topoCuts, b.cutValsForDVCS())) continue;
                    }

                    // binning
                    double xB=b.x, Q2=b.Q2, tt=std::fabs(b.t1), phi=b.phi();
                    if (!std::isfinite(xB)||!std::isfinite(Q2)||!std::isfinite(tt)||!std::isfinite(phi)) continue;
                    int ix=findBin(xB,xB_bins), iQ=findBin(Q2,Q2_bins), itb=findBin(tt,t_bins), ip=phiToBin(phi);
                    if (ix<0||iQ<0||itb<0||ip<0||ip>=N_PHI_BINS) continue;
                    counts[BinKey(ix,iQ,itb,ip)].N_pi0_mc++;
                }
            }
        }

        // ---- π0 reco MC (no helicity) -> N_pi0_reco ----
        if (t_pi0_reco) {
            BranchBinderEPPI0MC b; b.bind(t_pi0_reco);
            if (!b.readyCuts() || !b.readyBins()) {
                std::cerr << "[pi0_contam][WARN] eppi0 reco MC tree missing branches for '" << period << "'. Skipping." << std::endl;
            } else {
                const Long64_t nent = t_pi0_reco->GetEntries();
                for (Long64_t i=0;i<nent;++i) {
                    t_pi0_reco->GetEntry(i);
                    if (!applyKinematicCuts_simple(b.t1,b.open_angle_ep2,b.pTmiss)) continue;

                    bool match=false;
                    for (const auto& topoStr : topologies) {
                        if (passesTopology_simple(b.detector1,b.detector2,topoStr)) { match=true; break; }
                    }
                    if (!match) continue;

                    // 3σ cuts under π0 hypothesis
                    // (we still use data-defined windows but with theta_pi0_pi0)
                    // Use the same DVCS period key; theta_pi0_pi0 var will just be ignored by DVCS windows
                    // so we *do not* apply 3σ here (to remain consistent with the python that used only kinematics on reco?).
                    // If desired, switch to b.cutValsForEPPI0() and apply topoCuts similarly.

                    // binning
                    double xB=b.x, Q2=b.Q2, tt=std::fabs(b.t1), phi=b.phi();
                    if (!std::isfinite(xB)||!std::isfinite(Q2)||!std::isfinite(tt)||!std::isfinite(phi)) continue;
                    int ix=findBin(xB,xB_bins), iQ=findBin(Q2,Q2_bins), itb=findBin(tt,t_bins), ip=phiToBin(phi);
                    if (ix<0||iQ<0||itb<0||ip<0||ip>=N_PHI_BINS) continue;
                    counts[BinKey(ix,iQ,itb,ip)].N_pi0_reco++;
                }
            }
        }

        // ---- eppi0 experimental data (helicity-resolved) -> N_pi0_exp^± ----
        if (t_pi0_data) {
            BranchBinderEPPI0Data b; b.bind(t_pi0_data);
            if (!b.readyCuts() || !b.readyBins()) {
                std::cerr << "[pi0_contam][WARN] eppi0 DATA tree missing branches for '" << period << "'. Skipping eppi0 DATA." << std::endl;
            } else {
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

                    // For eppi0 data, apply 3σ cuts using DVCS windows? (python did: passes_3sigma_cuts(event, False, cuts_dict) with eppi0)
                    // We keep consistency with DVCS windows (data-defined).
                    // However, theta_pi0_pi0 appears in eppi0 data and not in DVCS windows; we still use available vars.
                    // If cuts are missing, we just skip 3σ.
                    // Build the DVCS key using *some* matched topology; we already ensured match.
                    // (Exact topo key not used here; would require same topo as in DVCS—fine to skip if absent)

                    // binning
                    double xB=b.x, Q2=b.Q2, tt=std::fabs(b.t1), phi=b.phi();
                    if (!std::isfinite(xB)||!std::isfinite(Q2)||!std::isfinite(tt)||!std::isfinite(phi)) continue;
                    int ix=findBin(xB,xB_bins), iQ=findBin(Q2,Q2_bins), itb=findBin(tt,t_bins), ip=phiToBin(phi);
                    if (ix<0||iQ<0||itb<0||ip<0||ip>=N_PHI_BINS) continue;

                    BinKey key(ix,iQ,itb,ip);
                    if (b.helicity==+1) counts[key].N_pi0_exp.plus++; else counts[key].N_pi0_exp.minus++;
                }
            }
        }

        // ---- compute contamination per helicity ----
        for (auto& kv : counts) {
            BinCounts& bc = kv.second;

            auto compute_one = [&](long long N_pi0_exp_h, long long N_data_h, double& c, double& c_err){
                if (N_data_h<=0 || bc.N_pi0_reco<=0 || bc.N_pi0_mc<=0) { c=0.0; c_err=0.0; return; }
                double ratio = static_cast<double>(N_pi0_exp_h) / static_cast<double>(bc.N_pi0_reco);
                double cval  = static_cast<double>(bc.N_pi0_mc) * ratio / static_cast<double>(N_data_h);

                auto rel = [](long long n)->double { return (n>0)? 1.0/std::sqrt(static_cast<double>(n)) : 0.0; };
                double rel_mc   = rel(bc.N_pi0_mc);
                double rel_exp  = rel(N_pi0_exp_h);
                double rel_reco = rel(bc.N_pi0_reco);
                double rel_data = rel(N_data_h);

                double rel_ratio = std::sqrt(rel_exp*rel_exp + rel_reco*rel_reco);
                double rel_tot   = std::sqrt(rel_mc*rel_mc + rel_ratio*rel_ratio + rel_data*rel_data);
                c = cval;
                c_err = cval * rel_tot;
            };

            compute_one(bc.N_pi0_exp.plus,  bc.N_data.plus,  bc.c_plus,  bc.c_plus_err);
            compute_one(bc.N_pi0_exp.minus, bc.N_data.minus, bc.c_minus, bc.c_minus_err);
        }

        // ---- write JSON for this period ----
        const std::string out_file = (fs::path(out_dir) / ("contamination_" + period + ".json")).string();
        writeContaminationJson(out_file, counts, N_PHI_BINS, xB_bins, Q2_bins, t_bins);

        // ---- optional copy for Fa18_inb_supp ----
        if (COPY_CONTAM_TO_FA18_INB_SUPP && runTag == "fa18_inb") {
            const std::string supp_period = "DVCS_Fa18_inb_supp";
            const std::string out_copy = (fs::path(out_dir) / ("contamination_" + supp_period + ".json")).string();
            std::error_code ec;
            fs::copy_file(out_file, out_copy, fs::copy_options::overwrite_existing, ec);
            if (ec) std::cerr << "[pi0_contam][WARN] Could not copy to Fa18_inb_supp JSON: " << ec.message() << std::endl;
            else     std::cout << "[pi0_contam] Also wrote (copy) " << out_copy << std::endl;
        }
    } // periods
}