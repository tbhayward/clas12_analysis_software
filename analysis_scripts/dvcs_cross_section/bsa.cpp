// bsa.cpp (uses pi0_corrected_counts_all_groups.json as the ONLY input for counts)
//
// Behavior summary (strict, with correct low-stat errors, no fallbacks):
// - Reads contamination-corrected helicity counts from:
//       <out_root_dir>/jsons/pi0_corrected_counts_all_groups.json
// - For each requested group in `periods`, computes per-bin BSA points and fits
//   the stabilized model: y(phi) = C + [A*sin(phi)] / [1 + B*cos(phi)]
// - Writes per-group fit JSONs to:  <out_root_dir>/jsons/BSA_fits/BSA_fits_<group>.json
// - Writes all-periods rollup to:   <out_root_dir>/jsons/BSA_fits_all_periods.json
// - Writes 10.6 combined to:        <out_root_dir>/jsons/BSA_fits_combined_10.6.json
// - Saves plots under:              <out_root_dir>/bsa_plots/<group>/plot_bsa_<group>_xB_<ix>.png
//
// Notes (names/conventions):
// - Corrected-master group keys are bare: "sp18_inb","sp18_out","fa18_inb","fa18_out","fa18_inb_supp","sp19_inb","10.6_GeV".
// - `periods` may contain either bare keys or DVCS_* keys (e.g. DVCS_Fa18_inb); we normalize internally.
// - Polarization P is now computed per (xB,Q2,t) cell ONLY (no phi subdivision), with NO extra exclusivity/topology
//   filters, then replicated to all phi bins. This avoids “counts exist but P missing” without touching upstream steps.
//
// Style: K&R braces, ASCII-only.

#include "bsa.h"

#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TH1.h>
#include <TF1.h>
#include <TLatex.h>
#include <TFitResult.h>
#include <TPad.h>
#include <TGaxis.h>
#include <TTree.h>

#include <Math/Factory.h>
#include <Math/Minimizer.h>
#include <Math/Functor.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace {

// =================== global style ===================
struct StyleInit {
    StyleInit() {
        gStyle->SetOptTitle(0);
        gStyle->SetOptStat(0);
        gStyle->SetFrameLineWidth(2);
        gStyle->SetLineWidth(2);
        gStyle->SetPadTickX(1);
        gStyle->SetPadTickY(1);
        gStyle->SetLegendBorderSize(1);
        const int rf = 42; // Helvetica
        gStyle->SetTitleFont(rf, "XYZ");
        gStyle->SetLabelFont(rf, "XYZ");
        gStyle->SetTextFont(rf);
    }
} _style_bootstrap;

constexpr int    N_PHI_BINS = 12;
constexpr double TWO_PI     = 2.0 * M_PI;

// ---- BSA fit stabilization knobs (single-cos denominator) ----
constexpr double EPS_DEN_FLOOR = 1e-2;   // demand D(phi) >= this across [0,2pi)
constexpr double LAMBDA_DEN    = 1e6;    // penalty multiplier for denominator floor
constexpr double EPS_DEN_EVAL  = 1e-6;   // evaluation clamp for denominator

constexpr double A_MAX_AMP     = 0.999;  // soft box for |B|
constexpr double LAMBDA_AMP    = 1e5;    // penalty multiplier for amplitude box

// Jeffreys prior alpha for counts
constexpr double JEFFREYS_ALPHA = 0.5;

using BinKey = std::tuple<int,int,int,int>; // (ix,iQ,it,ip)
struct HelCorr { double Np=0.0; double Nm=0.0; double ep=0.0; double em=0.0; };
using GroupTable = std::map<BinKey, HelCorr>;
using AllGroups  = std::map<std::string, GroupTable>;

// ------------ tiny helpers ------------
[[noreturn]] static void fatal(const std::string& msg) {
    std::cerr << "[bsa][FATAL] " << msg << std::endl;
    std::exit(EXIT_FAILURE);
}

static inline std::string key4s(int ix,int iQ,int it,int ip) {
    std::ostringstream os; os << "(" << ix << "," << iQ << "," << it << "," << ip << ")";
    return os.str();
}

static bool parse_tuple_key(const std::string& s, BinKey& out) {
    int ix,iQ,it,ip;
    if (std::sscanf(s.c_str(),"(%d,%d,%d,%d)",&ix,&iQ,&it,&ip)!=4) return false;
    out = BinKey(ix,iQ,it,ip);
    return true;
}

static inline std::vector<std::pair<double,double>> uniqueRanges(const std::vector<Binning>& scheme, char which) {
    std::set<std::pair<double,double>> s;
    for (const auto& b : scheme) {
        if (which == 'x') s.emplace(b.xBmin, b.xBmax);
        else if (which == 'Q') s.emplace(b.Q2min, b.Q2max);
        else if (which == 't') s.emplace(b.tmin, b.tmax);
    }
    return {s.begin(), s.end()};
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

static inline int findIndex(const std::pair<double,double>& range,
                            const std::vector<std::pair<double,double>>& ranges) {
    for (int i=0;i<(int)ranges.size();++i) if (ranges[i]==range) return i;
    return -1;
}

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

// Extract {...} block that follows a key (very small JSON helper)
static std::string objForKey(const std::string& s, const std::string& key) {
    size_t p = s.find(key);
    if (p==std::string::npos) fatal("Key '"+key+"' not found.");
    size_t br = s.find('{', p);
    if (br==std::string::npos) fatal("Malformed JSON after key '"+key+"'");
    int d=0; size_t i=br;
    for (; i<s.size(); ++i) {
        if (s[i]=='{') d++;
        else if (s[i]=='}') { d--; if (!d) { ++i; break; } }
    }
    if (d!=0) fatal("Unbalanced braces near key '"+key+"'");
    return s.substr(br, i-br);
}

static long long parseIntAfterColon_strict(const std::string& s, size_t cpos, const std::string& ctx) {
    size_t a = cpos + 1;
    while (a < s.size() && std::isspace((unsigned char)s[a])) ++a;
    size_t b = a;
    while (b < s.size() && (std::isdigit((unsigned char)s[b]) || s[b]=='-' || s[b]=='+')) ++b;
    try { return std::stoll(s.substr(a,b-a)); }
    catch(...) { fatal("Non-integer value in "+ctx); }
    return 0;
}

static double parseDoubleAfterColon_num(const std::string& s, size_t cpos, const std::string& ctx) {
    size_t a = cpos + 1;
    while (a < s.size() && std::isspace((unsigned char)s[a])) ++a;
    size_t b = a;
    auto isdelim = [](char c)->bool{
        return std::isspace((unsigned char)c) || c==',' || c=='}' || c==']';
    };
    while (b < s.size() && !isdelim(s[b])) ++b;
    std::string raw = s.substr(a, b-a);
    try { return std::stod(raw); } catch(...) { fatal("Non-numeric value in "+ctx); }
    return 0.0;
}

// ------------ load pi0_corrected_counts_all_groups.json ------------
struct BinningMeta { int phi_bins=0; size_t nx=0, nQ=0, nt=0; };

static BinningMeta load_meta_from_master(const std::string& s) {
    std::string meta = objForKey(s, "\"binning_meta\"");
    auto findN = [&](const char* k, const char* ctx)->int{
        size_t p = meta.find(k); if (p==std::string::npos) fatal(std::string("binning_meta missing ")+k);
        size_t c = meta.find(':', p); if (c==std::string::npos) fatal(std::string("binning_meta malformed for ")+k);
        return (int)parseIntAfterColon_strict(meta, c, ctx);
    };
    BinningMeta bm;
    bm.phi_bins = findN("\"phi_bins\"", "phi_bins");
    bm.nx       = (size_t)findN("\"xB_bins\"", "xB_bins");
    bm.nQ       = (size_t)findN("\"Q2_bins\"", "Q2_bins");
    bm.nt       = (size_t)findN("\"t_bins\"",  "t_bins");
    if (bm.phi_bins != N_PHI_BINS) fatal("phi_bins mismatch: expected 12.");
    return bm;
}

static AllGroups load_corrected_master(const std::string& master_path,
                                       BinningMeta& out_meta,
                                       std::vector<std::string>& group_order) {
    std::ifstream ifs(master_path);
    if (!ifs) fatal(std::string("Cannot open master corrected counts: ")+master_path);
    std::string s((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());

    out_meta = load_meta_from_master(s);

    std::string groupsObj = objForKey(s, "\"groups\"");

    AllGroups out;
    size_t gk=0;
    while (true) {
        size_t q1 = groupsObj.find('"', gk); if (q1==std::string::npos) break;
        size_t q2 = groupsObj.find('"', q1+1); if (q2==std::string::npos) fatal("Malformed group name.");
        std::string gname = groupsObj.substr(q1+1, q2-q1-1);
        group_order.push_back(gname);

        size_t binsS = groupsObj.find("\"bins\"", q2);
        if (binsS==std::string::npos) fatal("Group "+gname+" missing 'bins' object.");
        int d=0; size_t br = groupsObj.find('{', binsS); if (br==std::string::npos) fatal("Malformed 'bins' in group "+gname);
        size_t i=br; for (; i<groupsObj.size(); ++i) {
            if (groupsObj[i]=='{') ++d;
            else if (groupsObj[i]=='}') { --d; if (!d) { ++i; break; } }
        }
        std::string binsObj = groupsObj.substr(br, i-br);

        GroupTable tbl;
        size_t p=0; int nb=0;
        while (true) {
            size_t k1=binsObj.find('"', p); if (k1==std::string::npos) break;
            size_t k2=binsObj.find('"', k1+1); if (k2==std::string::npos) fatal("Malformed bin key in "+gname);
            std::string key = binsObj.substr(k1+1, k2-k1-1);
            BinKey bk; if (!parse_tuple_key(key, bk)) fatal("Bad bin tuple key in "+gname);

            size_t vS = binsObj.find('{', k2); if (vS==std::string::npos) fatal("Missing bin object in "+gname);
            int d2=0; size_t j=vS;
            for (; j<binsObj.size(); ++j) {
                if (binsObj[j]=='{') ++d2;
                else if (binsObj[j]=='}') { --d2; if (!d2) { ++j; break; } }
            }
            std::string v = binsObj.substr(vS, j-vS);

            size_t hp = v.find("\"+1\"");
            if (hp==std::string::npos) fatal("Missing +1 block in "+gname);
            size_t hv = v.find("\"value\"", hp); if (hv==std::string::npos) fatal("Missing +1.value in "+gname);
            size_t hc = v.find(':', hv); if (hc==std::string::npos) fatal("Malformed +1.value in "+gname);
            double Np = parseDoubleAfterColon_num(v, hc, gname+" (+1).value");

            size_t he = v.find("\"err\"", hp); if (he==std::string::npos) fatal("Missing +1.err in "+gname);
            size_t hce = v.find(':', he); if (hce==std::string::npos) fatal("Malformed +1.err in "+gname);
            double ep = parseDoubleAfterColon_num(v, hce, gname+" (+1).err");

            size_t mp = v.find("\"-1\"");
            if (mp==std::string::npos) fatal("Missing -1 block in "+gname);
            size_t mv = v.find("\"value\"", mp); if (mv==std::string::npos) fatal("Missing -1.value in "+gname);
            size_t mc = v.find(':', mv); if (mc==std::string::npos) fatal("Malformed -1.value in "+gname);
            double Nm = parseDoubleAfterColon_num(v, mc, gname+" (-1).value");

            size_t me = v.find("\"err\"", mp); if (me==std::string::npos) fatal("Missing -1.err in "+gname);
            size_t mce = v.find(':', me); if (mce==std::string::npos) fatal("Malformed -1.err in "+gname);
            double em = parseDoubleAfterColon_num(v, mce, gname+" (-1).err");

            (void)ep; (void)em;
            tbl[bk] = HelCorr{Np, Nm, ep, em};
            p = j; ++nb;
        }
        if (nb==0) {
            std::cerr << "[bsa][WARN] Group '"<<gname<<"' has zero bins in corrected master.\n";
        }
        out[gname] = std::move(tbl);
        gk = i;
    }
    if (out.empty()) fatal("No groups parsed from corrected master.");
    return out;
}

// ------------ name normalization / tree aliasing ------------
static const std::map<std::string,std::string> BARE_TO_DVCS = {
    {"sp18_inb",      "DVCS_Sp18_inb"},
    {"sp18_out",      "DVCS_Sp18_out"},
    {"fa18_inb_supp", "DVCS_Fa18_inb_supp"},
    {"fa18_inb",      "DVCS_Fa18_inb"},
    {"fa18_out",      "DVCS_Fa18_out"},
    {"sp19_inb",      "DVCS_Sp19_inb"}
};

static std::string to_bare_group(const std::string& name) {
    if (BARE_TO_DVCS.count(name)) return name;
    for (const auto& kv : BARE_TO_DVCS) if (kv.second == name) return kv.first;
    if (name == "10.6_GeV") return name;
    if (name.rfind("DVCS_", 0) == 0) {
        std::string tail = name.substr(5);
        std::transform(tail.begin(), tail.end(), tail.begin(),
                       [](unsigned char c){ return std::tolower(c); });
        return tail;
    }
    return name;
}

static TTree* find_tree_for_group(const std::map<std::string,TTree*>& trees,
                                  const std::string& requested_key,
                                  const std::string& bare_key) {
    auto it = trees.find(requested_key);
    if (it != trees.end() && it->second) return it->second;
    auto itAlias = BARE_TO_DVCS.find(bare_key);
    if (itAlias != BARE_TO_DVCS.end()) {
        auto it2 = trees.find(itAlias->second);
        if (it2 != trees.end() && it2->second) return it2->second;
    }
    auto itBare = trees.find(bare_key);
    if (itBare != trees.end() && itBare->second) return itBare->second;
    return nullptr;
}

// ------------ polarization (from DVCS tree), per (xB,Q2,t) cell ------------
struct PolStats { double P=0.0; int n=0; };

struct BranchPol {
    int helicity=0;
    double t1=0, open_angle_ep2=0, pTmiss=0;
    double x=0, Q2=0, phi2=0, Delta_phi=0;
    double beam_pol=1.0;
    int detector1=0, detector2=0;
    bool hasX=false, hasQ=false, hasT=false, hasPhi2=false, hasDp=false, hasPol=false;
    void bind(TTree* t){
        if (!t) return;
        auto bindI=[&](const char*n,int*a,bool&f){ if(t->GetBranch(n)){ t->SetBranchAddress(n,a); f=true; } };
        auto bindD=[&](const char*n,double*a,bool&f){ if(t->GetBranch(n)){ t->SetBranchAddress(n,a); f=true; } };
        bindD("x",&x,hasX);
        bindD("Q2",&Q2,hasQ);
        bindD("t1",&t1,hasT);
        bindD("phi2",&phi2,hasPhi2);
        bindD("Delta_phi",&Delta_phi,hasDp);
        bindD("beam_pol",&beam_pol,hasPol);
        (void)helicity; (void)open_angle_ep2; (void)pTmiss; (void)detector1; (void)detector2;
    }
    double phi() const { return hasPhi2 ? phi2 : (hasDp ? Delta_phi : std::numeric_limits<double>::quiet_NaN()); }
};

// Cell-only polarization (no phi subdivision, no extra cuts). Replicate to all phi later.
static std::vector<PolStats> compute_cell_polarization(
    TTree* t,
    const std::vector<std::pair<double,double>>& xB_bins,
    const std::vector<std::pair<double,double>>& Q2_bins,
    const std::vector<std::pair<double,double>>& t_bins)
{
    std::vector<PolStats> Pcell(xB_bins.size()*Q2_bins.size()*t_bins.size(), PolStats{0.0,0});
    auto idx3=[&](int ix,int iQ,int it){ return (ix*(int)Q2_bins.size()+iQ)*(int)t_bins.size()+it; };

    BranchPol b; b.bind(t);
    if (!(b.hasX && b.hasQ && b.hasT && b.hasPol)) return Pcell;

    const Long64_t nent = t->GetEntries();
    for (Long64_t i=0;i<nent;++i){
        t->GetEntry(i);
        double xB=b.x, Q2=b.Q2, tt=fabs(b.t1);
        if (!std::isfinite(xB)||!std::isfinite(Q2)||!std::isfinite(tt)) continue;

        auto findBin=[&](double v,const std::vector<std::pair<double,double>>& ranges)->int{
            for (int k=0;k<(int)ranges.size();++k) if (v>=ranges[k].first && v<ranges[k].second) return k;
            return -1;
        };
        int ix=findBin(xB,xB_bins), iQ=findBin(Q2,Q2_bins), it=findBin(tt,t_bins);
        if (ix<0||iQ<0||it<0) continue;

        PolStats& ps = Pcell[idx3(ix,iQ,it)];
        ps.P += (b.beam_pol - ps.P)/(ps.n+1);
        ps.n++;
    }
    return Pcell;
}

// ------------ math helpers for stabilizer penalties ------------
static inline double Dmin_single(double B) { return 1.0 - std::fabs(B); }
static inline double overBox(double A){ double v = std::fabs(A) - A_MAX_AMP; return (v>0.0)? v*v : 0.0; }

// ------------ BSA computation core ------------
struct BSApt { double phi=0.0; double bsa=0.0; double err=0.0; bool valid=false; };
struct FitRes { double A=0, Aerr=0, B1=0, B1err=0, B2=0, B2err=0, C=0, Cerr=0; double chi2=0; int ndf=0; int status=0; };
struct CellResult {
    std::vector<BSApt> points; FitRes fit;
    double P_used=std::numeric_limits<double>::quiet_NaN(); bool P_per_bin=true;
};

static FitRes fit_cell(const std::vector<BSApt>& pts){
    std::vector<double> x, y, ey;
    x.reserve(pts.size()); y.reserve(pts.size()); ey.reserve(pts.size());
    for (const auto& p : pts) {
        if (!p.valid) continue;
        x.push_back(p.phi);
        y.push_back(p.bsa);
        ey.push_back(std::max(p.err, 1e-6));
    }
    FitRes fr; 
    const int n = (int)x.size();
    if (n < 4) { fr.status = 1; fr.ndf = 0; return fr; }

    auto chi2pen = [&](const double *par){
        const double A  = par[0];
        const double B1 = par[1];
        const double C  = par[2];

        double chi2 = 0.0;
        for (int i=0;i<n;++i){
            const double phi = x[i];
            double denom = 1.0 + B1*std::cos(phi);
            if (denom < EPS_DEN_EVAL) denom = EPS_DEN_EVAL;
            const double yhat = C + (A*std::sin(phi))/denom;
            const double pull = (y[i] - yhat)/ey[i];
            chi2 += pull*pull;
        }

        const double Dmin = Dmin_single(B1);
        double pen_den = 0.0;
        if (Dmin < EPS_DEN_FLOOR) {
            const double deficit = EPS_DEN_FLOOR - Dmin;
            pen_den = LAMBDA_DEN * deficit * deficit;
        }
        double pen_amp = overBox(B1);
        return chi2 + pen_den + LAMBDA_AMP*pen_amp;
    };

    std::unique_ptr<ROOT::Math::Minimizer> min(
        ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"));
    min->SetMaxFunctionCalls(10000);
    min->SetMaxIterations(10000);
    min->SetTolerance(1e-6);
    min->SetPrintLevel(0);

    ROOT::Math::Functor fcn(chi2pen, 3);
    min->SetFunction(fcn);

    min->SetLimitedVariable(0, "A",  0.10, 0.02, -1.0,  1.0);
    min->SetLimitedVariable(1, "B1", 0.00, 0.02, -1.0,  1.0);
    min->SetLimitedVariable(2, "C",  0.00, 0.02, -0.5,  0.5);

    const bool ok = min->Minimize();
    fr.status = ok ? 0 : 1;

    const double *par = min->X();
    const double *err = min->Errors();

    fr.A  = par[0]; fr.Aerr  = err[0];
    fr.B1 = par[1]; fr.B1err = err[1];
    fr.B2 = 0.0;    fr.B2err = 0.0;
    fr.C  = par[2]; fr.Cerr  = err[2];

    fr.chi2 = chi2pen(par);
    fr.ndf  = std::max(0, n - 3);
    return fr;
}

// JSON writers (schema preserved)
static void write_period_bsa_json(
    const std::string& out_path,
    int nPhi,
    const std::vector<std::pair<double,double>>& xB_bins,
    const std::vector<std::pair<double,double>>& Q2_bins,
    const std::vector<std::pair<double,double>>& t_bins,
    const std::map<std::tuple<int,int,int>, CellResult>& cells)
{
    std::ofstream ofs(out_path);
    if (!ofs) { std::cerr<<"[bsa][ERROR] Cannot open "<<out_path<<"\n"; return; }
    ofs<<std::fixed<<std::setprecision(8);
    ofs<<"{\n";
    ofs<<"  \"binning_meta\": {\"phi_bins\": "<<nPhi<<", \"xB_bins\": "<<xB_bins.size()<<", \"Q2_bins\": "<<Q2_bins.size()<<", \"t_bins\": "<<t_bins.size()<<"},\n";
    ofs<<"  \"bins\": {\n";
    bool first=true;
    for (const auto& kv : cells){
        if (!first) ofs<<",\n"; first=false;
        int ix,iQ,it; std::tie(ix,iQ,it)=kv.first;
        const auto& cr = kv.second;
        ofs<<"    \"("<<ix<<","<<iQ<<","<<it<<")\": {\n";
        ofs<<"      \"phi\": [";
        for (size_t i=0;i<cr.points.size();++i){ if(i) ofs<<","; ofs<<cr.points[i].phi; } ofs<<"],\n";
        ofs<<"      \"bsa\": [";
        for (size_t i=0;i<cr.points.size();++i){ if(i) ofs<<","; ofs<<cr.points[i].bsa; } ofs<<"],\n";
        ofs<<"      \"bsa_err\": [";
        for (size_t i=0;i<cr.points.size();++i){ if(i) ofs<<","; ofs<<cr.points[i].err; } ofs<<"],\n";
        ofs<<"      \"fit\": {"
              "\"A\": {\"val\": "<<cr.fit.A<<", \"err\": "<<cr.fit.Aerr<<"}, "
              "\"B1\": {\"val\": "<<cr.fit.B1<<", \"err\": "<<cr.fit.B1err<<"}, "
              "\"B2\": {\"val\": 0.0, \"err\": 0.0}, "
              "\"C\": {\"val\": "<<cr.fit.C<<", \"err\": "<<cr.fit.Cerr<<"}, "
              "\"chi2\": "<<cr.fit.chi2<<", \"ndf\": "<<cr.fit.ndf<<", \"status\": "<<cr.fit.status<<"},\n";
        ofs<<"      \"polarization\": {\"per_bin\": "<<(cr.P_per_bin?"true":"false")<<", \"P_used\": "<<(std::isfinite(cr.P_used)? cr.P_used : 0.0)<<"}\n";
        ofs<<"    }";
    }
    ofs<<"\n  }\n}\n";
    std::cout << "[bsa] Wrote " << out_path << "\n";
}

static void write_all_periods_json(
    const std::string& out_path,
    const std::map<std::string, std::map<std::tuple<int,int,int>, CellResult>>& perPeriod,
    int nPhi,
    const std::vector<std::pair<double,double>>& xB_bins,
    const std::vector<std::pair<double,double>>& Q2_bins,
    const std::vector<std::pair<double,double>>& t_bins)
{
    std::ofstream ofs(out_path);
    if (!ofs) { std::cerr<<"[bsa][ERROR] Cannot open "<<out_path<<"\n"; return; }
    ofs<<std::fixed<<std::setprecision(8);
    ofs<<"{\n";
    ofs<<"  \"binning_meta\": {\"phi_bins\": "<<nPhi<<", \"xB_bins\": "<<xB_bins.size()<<", \"Q2_bins\": "<<Q2_bins.size()<<", \"t_bins\": "<<t_bins.size()<<"},\n";
    ofs<<"  \"periods\": {\n";
    bool firstP=true;
    for (const auto& pkv : perPeriod){
        if (!firstP) ofs<<",\n"; firstP=false;
        ofs<<"    \""<<pkv.first<<"\": {\n      \"bins\": {\n";
        bool firstB=true;
        for (const auto& kv : pkv.second){
            if (!firstB) ofs<<",\n"; firstB=false;
            int ix,iQ,it; std::tie(ix,iQ,it)=kv.first;
            const auto& cr = kv.second;
            ofs<<"        \"("<<ix<<","<<iQ<<","<<it<<")\": {";
            ofs<<"\"phi\":["; for (size_t i=0;i<cr.points.size();++i){ if(i)ofs<<","; ofs<<cr.points[i].phi; } ofs<<"],";
            ofs<<"\"bsa\":["; for (size_t i=0;i<cr.points.size();++i){ if(i)ofs<<","; ofs<<cr.points[i].bsa; } ofs<<"],";
            ofs<<"\"bsa_err\":["; for (size_t i=0;i<cr.points.size();++i){ if(i)ofs<<","; ofs<<cr.points[i].err; } ofs<<"],";
            ofs<<"\"fit\":{"
                  "\"A\":{\"val\":"<<cr.fit.A<<",\"err\":"<<cr.fit.Aerr<<"},"
                  "\"B1\":{\"val\":"<<cr.fit.B1<<",\"err\":"<<cr.fit.B1err<<"},"
                  "\"B2\":{\"val\":0.0,\"err\":0.0},"
                  "\"C\":{\"val\":"<<cr.fit.C<<",\"err\":"<<cr.fit.Cerr<<"},"
                  "\"chi2\":"<<cr.fit.chi2<<",\"ndf\":"<<cr.fit.ndf<<",\"status\":"<<cr.fit.status<<"},";
            ofs<<"\"polarization\":{\"per_bin\":"<<(cr.P_per_bin?"true":"false")<<",\"P_used\":"<<(std::isfinite(cr.P_used)? cr.P_used : 0.0)<<"}";
            ofs<<"}";
        }
        ofs<<"\n      }\n    }";
    }
    ofs<<"\n  }\n}\n";
    std::cout << "[bsa] Wrote " << out_path << "\n";
}

// ------------ plotting ------------
static void drawDegreeTicks(double xmin, double ymin, double xmax, double labelSize){
    TGaxis* ax = new TGaxis(xmin, ymin, xmax, ymin, 0.0, 360.0, 4, "");
    ax->SetLabelFont(42);
    ax->SetLabelSize(labelSize);
    ax->SetLabelOffset(0.012);
    ax->SetTitle("");
    ax->SetTickSize(0.02);
    ax->Draw();
}

static std::string plot_subdir_for_group(const std::string& group) {
    if (group == "10.6_GeV") return "10.6_combined";
    return group;
}

static void plot_cells_for_period(
    const std::string& period,
    const std::vector<Binning>& binning_scheme,
    const std::vector<std::pair<double,double>>& xB_bins,
    const std::vector<std::pair<double,double>>& Q2_bins,
    const std::vector<std::pair<double,double>>& t_bins,
    const std::map<std::tuple<int,int,int>, CellResult>& cells,
    const std::string& out_dir_plots)
{
    namespace fs = std::filesystem;
    std::error_code ec;
    fs::create_directories(out_dir_plots, ec);

    static const auto PHI_DEG = phiCentersDeg();

    for (int ix = 0; ix < (int)xB_bins.size(); ++ix) {
        const auto xb = xB_bins[ix];

        std::vector<std::pair<double,double>> Q2_slice, t_slice;
        uniqueQT_for_xB(binning_scheme, xb, Q2_slice, t_slice);
        if (Q2_slice.empty() || t_slice.empty()) continue;

        const int nrows = (int)t_slice.size();
        const int ncols = (int)Q2_slice.size();

        const int W = 280*ncols + 160;
        const int H = 240*nrows + 170;

        std::ostringstream cname; cname<<"c_bsa_"<<period<<"_xB"<<ix;
        TCanvas* c = new TCanvas(cname.str().c_str(), cname.str().c_str(), W, H);

        TPad* pTop  = new TPad("pTop","pTop", 0.0, 0.915, 1.0, 1.0);
        pTop->SetFillStyle(0); pTop->SetBorderSize(0); pTop->Draw();

        TPad* pGrid = new TPad("pGrid","pGrid", 0.0, 0.00, 1.0, 0.915);
        pGrid->SetFillStyle(0); pGrid->SetBorderSize(0); pGrid->Draw();
        pGrid->cd();
        pGrid->Divide(ncols, nrows, 0.0001, 0.0001);

        pTop->cd();
        TLatex head;
        head.SetNDC(); head.SetTextAlign(22);
        head.SetTextFont(42);
        head.SetTextSize(0.36);
        std::ostringstream tit;
        tit << Form("Beam-Spin Asymmetry  %s   x_{B} #in [%.2g, %.2g]",
                    period.c_str(), xb.first, xb.second);
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
                gPad->SetLeftMargin(0.125);
                gPad->SetRightMargin(0.10);

                TH1* frame = gPad->DrawFrame(0.0, -1.05, 360.0, 1.05);
                TAxis* ax = frame->GetXaxis();
                TAxis* ay = frame->GetYaxis();

                ax->SetLabelSize(0.0001);
                ax->SetTitle("#phi (deg)");
                ay->SetTitle("A_{LU}");
                ax->CenterTitle(); ay->CenterTitle();
                ax->SetNdivisions(505);

                ax->SetTitleSize(0.060);
                ay->SetTitleSize(0.060);
                ay->SetLabelSize(0.048);

                ax->SetTitleOffset(1.25);
                ay->SetTitleOffset(1.35);

                drawDegreeTicks(gPad->GetUxmin(), gPad->GetUymin(), gPad->GetUxmax(), 0.048);

                auto itCell = cells.find(std::make_tuple(ix, iQ_global, it_global));
                if (itCell == cells.end()) continue;
                const auto& cr = itCell->second;

                std::vector<double> x, y, ey;
                x.reserve(N_PHI_BINS); y.reserve(N_PHI_BINS); ey.reserve(N_PHI_BINS);
                for (int ip=0; ip<N_PHI_BINS; ++ip) {
                    const auto& p = cr.points[ip];
                    if (!p.valid) continue;
                    x.push_back(PHI_DEG[ip]);
                    y.push_back(p.bsa);
                    ey.push_back(std::max(1e-6, p.err));
                }
                if (!x.empty()) {
                    TGraphErrors* gr = new TGraphErrors((int)x.size(), x.data(), y.data(), nullptr, ey.data());
                    gr->SetMarkerStyle(20);
                    gr->SetMarkerSize(1.1);
                    gr->SetLineWidth(2);
                    gr->Draw("P SAME");
                }

                if (cr.fit.status == 0 || cr.fit.ndf > 0) {
                    const int NS=721;
                    std::vector<double> xd(NS), yd(NS);
                    for (int i=0;i<NS;++i){
                        double deg = double(i)*0.5;
                        double rad = deg * (TWO_PI/360.0);
                        double denom = 1.0 + cr.fit.B1*std::cos(rad);
                        if (denom < EPS_DEN_EVAL) denom = EPS_DEN_EVAL;
                        double val = cr.fit.C + (cr.fit.A*std::sin(rad))/denom;
                        xd[i] = deg; yd[i] = val;
                    }
                    TGraph* gfit = new TGraph(NS, xd.data(), yd.data());
                    gfit->SetLineColor(kRed);
                    gfit->SetLineWidth(2);
                    gfit->Draw("L SAME");
                }

                TLatex lab;
                lab.SetNDC(); lab.SetTextSize(0.040); lab.SetTextAlign(11);
                lab.SetTextFont(42);
                lab.DrawLatex(0.15, 0.94,
                    Form("Q^{2} #in [%.2g, %.2g],   -t #in [%.2g, %.2g]",
                         Q2_slice[ccol].first, Q2_slice[ccol].second,
                         t_slice[r].first,    t_slice[r].second));

                TLegend* leg = new TLegend(0.50, 0.68, 0.90, 0.92);
                leg->SetBorderSize(1);
                leg->SetLineColor(kBlack);
                leg->SetFillColor(kWhite);
                leg->SetFillStyle(1001);
                leg->SetTextFont(42);
                leg->SetTextSize(0.040);
                leg->AddEntry((TObject*)nullptr, Form("A = %.3f +/- %.3f",  cr.fit.A,  cr.fit.Aerr), "");
                leg->AddEntry((TObject*)nullptr, Form("B = %.3f +/- %.3f",  cr.fit.B1, cr.fit.B1err), "");
                leg->AddEntry((TObject*)nullptr, Form("C = %.3f +/- %.3f",  cr.fit.C,  cr.fit.Cerr), "");
                leg->Draw();
            }
        }

        std::ostringstream fout;
        fout << out_dir_plots << "/plot_bsa_" << period << "_xB_" << ix << ".png";
        c->SaveAs(fout.str().c_str());
        delete c;
        std::cout << "[bsa] Wrote " << fout.str() << "\n";
    }
}

} // end anonymous namespace

// =======================================================
// Public driver
// =======================================================
void compute_and_plot_bsa_helicity(
    const std::vector<std::string>& periods,
    const std::vector<std::string>& topologies, // unused in P computation now (kept for signature stability)
    const std::vector<Binning>& binning_scheme,
    const std::map<std::string, TTree*>& dvcsDataTrees,
    const std::string& pi0_corrected_counts_json_path,
    const std::string& out_root_dir)
{
    namespace fs = std::filesystem;

    // Load corrected counts (all groups, bare keys)
    BinningMeta master_meta;
    std::vector<std::string> group_order_in_master;
    AllGroups allGroups = load_corrected_master(pi0_corrected_counts_json_path, master_meta, group_order_in_master);

    // Build binning from your runtime scheme (checks sizes)
    const auto xB_bins = uniqueRanges(binning_scheme, 'x');
    const auto Q2_bins = uniqueRanges(binning_scheme, 'Q');
    const auto t_bins  = uniqueRanges(binning_scheme, 't');
    if (xB_bins.size() != master_meta.nx || Q2_bins.size() != master_meta.nQ || t_bins.size() != master_meta.nt) {
        fatal("Binning mismatch between runtime binning_scheme and corrected master binning_meta.");
    }

    const auto PHI_RAD = phiCentersRad();

    const fs::path json_period_dir = fs::path(out_root_dir)/"jsons"/"BSA_fits";
    std::error_code ecDir;
    fs::create_directories(json_period_dir, ecDir);

    std::map<std::string, std::map<std::tuple<int,int,int>, CellResult>> allPeriodCells;

    for (const auto& req_name : periods) {
        const std::string bare = to_bare_group(req_name);

        auto itG = allGroups.find(bare);
        if (itG == allGroups.end()) {
            fatal("Requested group '"+req_name+"' (normalized to '"+bare+"') not present in corrected master.");
        }
        const GroupTable& table = itG->second;

        TTree* t = find_tree_for_group(dvcsDataTrees, req_name, bare);
        if (!t) {
            std::ostringstream os;
            os << "Missing DVCS polarization tree for requested group '"<<req_name<<"' "
               << "(normalized to '"<<bare<<"'). Tried exact key, canonical DVCS alias, and bare key.";
            fatal(os.str());
        }

        // Compute polarization per (xB,Q2,t) cell; replicate to all phi bins.
        std::vector<PolStats> Pcell = compute_cell_polarization(t, xB_bins, Q2_bins, t_bins);
        auto idx3=[&](int ix,int iQ,int it){ return (ix*(int)Q2_bins.size()+iQ)*(int)t_bins.size()+it; };

        std::map<std::tuple<int,int,int>, CellResult> cells;
        for (int ix=0; ix<(int)xB_bins.size(); ++ix)
        for (int iQ=0; iQ<(int)Q2_bins.size(); ++iQ)
        for (int itb=0; itb<(int)t_bins.size(); ++itb) {
            CellResult result;
            result.points.resize(N_PHI_BINS);
            result.P_per_bin = false; // now per-cell
            const size_t f3 = idx3(ix,iQ,itb);
            const double P_here = (f3 < Pcell.size() && Pcell[f3].n > 0) ? Pcell[f3].P : 0.0;
            result.P_used = P_here;

            for (int ip=0; ip<N_PHI_BINS; ++ip) {
                const BinKey bk(ix,iQ,itb,ip);

                double Np = 0.0, Nm = 0.0;
                auto itB = table.find(bk);
                if (itB != table.end()) {
                    Np = itB->second.Np;
                    Nm = itB->second.Nm;
                }

                BSApt p; p.phi = PHI_RAD[ip];

                if (!(Np > 0.0 && Nm > 0.0)) {
                    p.valid = false; p.bsa=0.0; p.err=0.0;
                    result.points[ip] = p;
                    continue;
                }

                if (!(P_here > 0.0)) {
                    fatal("Polarization missing for group '"+bare+"' cell ("+
                          std::to_string(ix)+","+std::to_string(iQ)+","+std::to_string(itb)+
                          ") while counts exist (S>0).");
                }

                const double a = Np + JEFFREYS_ALPHA;
                const double b = Nm + JEFFREYS_ALPHA;
                const double S = a + b;
                const double D = a - b;

                const double A_raw   = D / S;
                const double var_raw = 4.0 * (a*b) / ((a+b)*(a+b)*(a+b+1.0));
                p.bsa  = A_raw / P_here;
                p.err  = std::sqrt(std::max(var_raw, 1e-12)) / P_here;
                p.valid = std::isfinite(p.bsa) && std::isfinite(p.err);

                result.points[ip] = p;
            }

            result.fit = fit_cell(result.points);
            cells[std::make_tuple(ix,iQ,itb)] = std::move(result);
        }

        const fs::path outP = json_period_dir/("BSA_fits_"+bare+".json");
        write_period_bsa_json(outP.string(), N_PHI_BINS, xB_bins, Q2_bins, t_bins, cells);

        const fs::path plots_dir = fs::path(out_root_dir)/"bsa_plots"/plot_subdir_for_group(bare);
        std::error_code ec; fs::create_directories(plots_dir, ec);
        plot_cells_for_period(bare, binning_scheme, xB_bins, Q2_bins, t_bins, cells, plots_dir.string());

        allPeriodCells[bare] = std::move(cells);
    }

    // all-periods rollup
    write_all_periods_json(
        (fs::path(out_root_dir)/"jsons"/"BSA_fits_all_periods.json").string(),
        allPeriodCells, N_PHI_BINS, xB_bins, Q2_bins, t_bins);

    // ---- Combined 10.6 output ----
    auto it106 = allGroups.find("10.6_GeV");
    if (it106 == allGroups.end()) fatal("No '10.6_GeV' group in corrected master; cannot build combined output.");

    const std::vector<std::string> comps_bare = {"sp18_inb","sp18_out","fa18_inb","fa18_out"};
    std::map<std::string,TTree*> compTrees;
    for (const auto& cg : comps_bare) {
        if (!allGroups.count(cg)) fatal("10.6_GeV requires component counts group '"+cg+"'.");
        TTree* tcomp = find_tree_for_group(dvcsDataTrees, cg, cg);
        if (!tcomp) fatal(std::string("Missing DVCS polarization tree for component '")+cg+"' for 10.6_GeV.");
        compTrees[cg] = tcomp;
    }

    // Per-component P per (xB,Q2,t), then weighted by component corrected counts
    std::map<std::string, std::vector<PolStats>> P_comp;
    for (const auto& kv : compTrees) {
        P_comp[kv.first] = compute_cell_polarization(kv.second, xB_bins, Q2_bins, t_bins);
    }

    std::map<std::tuple<int,int,int>, CellResult> combCells;
    const GroupTable& table106 = it106->second;
    auto idx3=[&](int ix,int iQ,int it){ return (ix*(int)Q2_bins.size()+iQ)*(int)t_bins.size()+it; };

    for (int ix=0; ix<(int)xB_bins.size(); ++ix)
    for (int iQ=0; iQ<(int)Q2_bins.size(); ++iQ)
    for (int itb=0; itb<(int)t_bins.size(); ++itb) {
        CellResult cr; cr.points.resize(N_PHI_BINS); cr.P_per_bin=false;

        // counts-weighted P across components (per cell; same for all phi)
        double wsum = 0.0, psum = 0.0;
        for (const auto& cg : comps_bare) {
            const auto& tblC = allGroups.at(cg);
            double S_cell = 0.0;
            for (int ip=0; ip<N_PHI_BINS; ++ip) {
                auto itC = tblC.find(BinKey(ix,iQ,itb,ip));
                if (itC != tblC.end() && itC->second.Np>0.0 && itC->second.Nm>0.0)
                    S_cell += itC->second.Np + itC->second.Nm;
            }
            const auto& Pc = P_comp[cg];
            const size_t f3 = idx3(ix,iQ,itb);
            const double Pk = (f3 < Pc.size() && Pc[f3].n > 0) ? Pc[f3].P : 0.0;
            if (S_cell > 0.0 && Pk > 0.0) { wsum += S_cell; psum += S_cell * Pk; }
        }
        if (!(wsum > 0.0)) {
            fatal("Could not build polarization for 10.6_GeV at cell ("+
                  std::to_string(ix)+","+std::to_string(iQ)+","+std::to_string(itb)+").");
        }
        const double P_here = psum/wsum;
        cr.P_used = P_here;

        for (int ip=0; ip<N_PHI_BINS; ++ip){
            const BinKey bk(ix,iQ,itb,ip);
            BSApt p; p.phi = phiCentersRad()[ip]; p.valid=false;

            double Np = 0.0, Nm = 0.0;
            auto it = table106.find(bk);
            if (it != table106.end()) { Np = it->second.Np; Nm = it->second.Nm; }

            if (!(Np > 0.0 && Nm > 0.0)) { cr.points[ip] = p; continue; }

            const double a = Np + JEFFREYS_ALPHA;
            const double b = Nm + JEFFREYS_ALPHA;
            const double S = a + b;
            const double D = a - b;

            const double A_raw   = D / S;
            const double var_raw = 4.0 * (a*b) / ((a+b)*(a+b)*(a+b+1.0));
            p.bsa  = A_raw / P_here;
            p.err  = std::sqrt(std::max(var_raw, 1e-12)) / P_here;
            p.valid = std::isfinite(p.bsa) && std::isfinite(p.err);

            cr.points[ip] = p;
        }

        cr.fit = fit_cell(cr.points);
        combCells[std::make_tuple(ix,iQ,itb)] = std::move(cr);
    }

    const fs::path outC = fs::path(out_root_dir)/"jsons"/"BSA_fits_combined_10.6.json";
    write_period_bsa_json(outC.string(), N_PHI_BINS, xB_bins, Q2_bins, t_bins, combCells);

    const fs::path plots_comb106 = fs::path(out_root_dir)/"bsa_plots"/"10.6_combined";
    std::error_code ec; fs::create_directories(plots_comb106, ec);
    plot_cells_for_period("RGA_10.6_combined", binning_scheme, xB_bins, Q2_bins, t_bins, combCells, plots_comb106.string());
}