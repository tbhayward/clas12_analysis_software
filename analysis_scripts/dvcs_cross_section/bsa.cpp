// bsa.cpp
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

// ---- BSA fit stabilization knobs (π+ style, adapted to BSA) ----
constexpr double EPS_DEN_FLOOR = 1e-2;   // demand D(φ) ≥ this over [0,2π)
constexpr double LAMBDA_DEN    = 1e6;    // penalty multiplier for denominator floor
constexpr double EPS_DEN_EVAL  = 1e-6;   // evaluation clamp for denominator

constexpr double A_MAX_AMP     = 0.999;  // soft "box" for |B1|,|B2| (and optional |A|,|C|)
constexpr double LAMBDA_AMP    = 1e5;    // penalty multiplier for amplitude box

using BinKey = std::tuple<int,int,int,int>; // (ix,iQ2,it,ip)
struct HelCounts { long long plus=0, minus=0; };

static inline std::string toLower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c){ return std::tolower(c); });
    return s;
}
static std::string periodToRunTagKey(const std::string& period) {
    auto pos = period.find('_');
    if (pos == std::string::npos || pos + 1 >= period.size()) return toLower(period);
    return toLower(period.substr(pos + 1));
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

// For a given xB range, list only the Q² and |t| ranges present in the CSV/binning.
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

// ------------ I/O helpers: total_counts.json ------------
using GroupCounts = std::map<std::string, std::map<BinKey, HelCounts>>;

static bool parse_tuple_key(const std::string& s, BinKey& out) {
    int ix,iQ,it,ip;
    if (std::sscanf(s.c_str(),"(%d,%d,%d,%d)",&ix,&iQ,&it,&ip)!=4) return false;
    out = BinKey(ix,iQ,it,ip);
    return true;
}

// Minimal parser (shape written by your total_counts code)
static bool load_total_counts(const std::string& path, GroupCounts& outGroups) {
    std::ifstream ifs(path);
    if (!ifs) { std::cerr<<"[bsa][ERROR] Cannot open total_counts JSON: "<<path<<"\n"; return false; }
    std::string s((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());

    size_t gpos = s.find("\"groups\"");
    if (gpos==std::string::npos) { std::cerr<<"[bsa][ERROR] 'groups' not found in total_counts.\n"; return false; }
    size_t brace = s.find('{', gpos); if (brace==std::string::npos) return false;
    int d=0; size_t i=brace; for (; i<s.size(); ++i){ if(s[i]=='{') d++; else if(s[i]=='}'){ d--; if(!d){ ++i; break; } } }
    std::string groupsObj = s.substr(brace, i-brace);

    size_t kpos=0;
    while (true) {
        size_t q1 = groupsObj.find('"', kpos); if (q1==std::string::npos) break;
        size_t q2 = groupsObj.find('"', q1+1); if (q2==std::string::npos) break;
        std::string gname = groupsObj.substr(q1+1, q2-q1-1);

        size_t objS = groupsObj.find('{', q2); if (objS==std::string::npos) break;
        int d2=0; size_t j=objS;
        for (; j<groupsObj.size(); ++j){ if(groupsObj[j]=='{') d2++; else if(groupsObj[j]=='}'){ d2--; if(!d2){ ++j; break; } } }
        std::string binsObj = groupsObj.substr(objS, j-objS);

        std::map<BinKey, HelCounts> gmap;
        size_t bpos=0;
        while (true) {
            size_t bk1 = binsObj.find('"', bpos); if (bk1==std::string::npos) break;
            size_t bk2 = binsObj.find('"', bk1+1); if (bk2==std::string::npos) break;
            std::string key = binsObj.substr(bk1+1, bk2-bk1-1);
            BinKey bk;
            if (!parse_tuple_key(key, bk)) { bpos=bk2+1; continue; }

            size_t valS = binsObj.find('{', bk2); if (valS==std::string::npos) break;
            int d3=0; size_t jj=valS;
            for (; jj<binsObj.size(); ++jj){ if(binsObj[jj]=='{') d3++; else if(binsObj[jj]=='}'){ d3--; if(!d3){ ++jj; break; } } }
            std::string obj = binsObj.substr(valS, jj-valS);

            auto findLL = [&](const char* pat)->long long {
                size_t p = obj.find(pat); if (p==std::string::npos) return 0;
                p = obj.find(':', p); if (p==std::string::npos) return 0;
                size_t a=p+1; while (a<obj.size() && isspace((unsigned char)obj[a])) ++a;
                size_t b=a; while (b<obj.size() && (isdigit((unsigned char)obj[b])||obj[b]=='-')) ++b;
                try { return std::stoll(obj.substr(a,b-a)); } catch(...) { return 0; }
            };
            HelCounts hc;
            hc.plus  = findLL("\"+1\"");
            hc.minus = findLL("\"-1\"");

            gmap[bk]=hc;
            bpos=jj;
        }

        outGroups[gname]=std::move(gmap);
        kpos=j;
    }
    return !outGroups.empty();
}

// ------------ I/O helpers: contamination_<period>.json ------------
struct ContamBin { double c_plus=0, c_plus_err=0, c_minus=0, c_minus_err=0; };
using ContamMap = std::map<BinKey, ContamBin>;

static bool load_contam_for_period(const std::string& path, ContamMap& out) {
    std::ifstream ifs(path);
    if (!ifs) return false;
    std::string s((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    size_t pos = s.find("\"bins\""); if (pos==std::string::npos) return false;
    size_t br = s.find('{', pos); if (br==std::string::npos) return false;
    int d=0; size_t i=br; for (; i<s.size(); ++i){ if(s[i]=='{') d++; else if(s[i]=='}'){ d--; if(!d){ ++i; break; } } }
    std::string binsObj = s.substr(br, i-br);

    size_t kpos=0;
    while (true) {
        size_t q1 = binsObj.find('"', kpos); if (q1==std::string::npos) break;
        size_t q2 = binsObj.find('"', q1+1); if (q2==std::string::npos) break;
        std::string key = binsObj.substr(q1+1, q2-q1-1);
        BinKey bk; if (!parse_tuple_key(key, bk)) { kpos=q2+1; continue; }

        size_t objS = binsObj.find('{', q2); if (objS==std::string::npos) break;
        int d2=0; size_t j=objS; for (; j<binsObj.size(); ++j){ if(binsObj[j]=='{') d2++; else if(binsObj[j]=='}'){ d2--; if(!d2){ ++j; break; } } }
        std::string obj = binsObj.substr(objS, j-objS);

        auto findD = [&](const char* pat)->double{
            size_t p=obj.find(pat); if (p==std::string::npos) return 0.0;
            p=obj.find(':',p); if (p==std::string::npos) return 0.0;
            size_t a=p+1; while (a<obj.size() && isspace((unsigned char)obj[a])) ++a;
            size_t b=a; while (b<obj.size() && (isdigit((unsigned char)obj[b])||obj[b]=='-'||obj[b]=='.'||obj[b]=='e'||obj[b]=='E'||obj[b]=='+')) ++b;
            try { return std::stod(obj.substr(a,b-a)); } catch(...) { return 0.0; }
        };
        ContamBin cb;
        cb.c_plus      = findD("\"contamination\":{\"+1\":{\"value\"");
        cb.c_plus_err  = findD("\"contamination\":{\"+1\":{\"err\"");
        cb.c_minus     = findD("\"contamination\":{\"-1\":{\"value\"");
        cb.c_minus_err = findD("\"contamination\":{\"-1\":{\"err\"");
        out[bk]=cb;
        kpos=j;
    }
    return true;
}

// ------------ polarization (from DVCS tree) ------------
struct PolStats { double P=1.0; double P_se=0.0; int n=0; };

struct BranchPol {
    int helicity=0;
    double t1=0, open_angle_ep2=0, pTmiss=0;
    double x=0, Q2=0, phi2=0, Delta_phi=0;
    double beam_pol=1.0;
    int detector1=0, detector2=0;
    bool hasH=false, hasCuts=false, hasBins=false, hasPhi2=false, hasDp=false, hasPol=false, hasTopo=false;
    void bind(TTree* t){
        if (!t) return;
        auto bindI=[&](const char*n,int*a,bool&f){ if(t->GetBranch(n)){ t->SetBranchAddress(n,a); f=true; } };
        auto bindD=[&](const char*n,double*a,bool&f){ if(t->GetBranch(n)){ t->SetBranchAddress(n,a); f=true; } };
        bindI("helicity",&helicity,hasH);
        bindD("t1",&t1,hasCuts); bindD("open_angle_ep2",&open_angle_ep2,hasCuts); bindD("pTmiss",&pTmiss,hasCuts);
        bindD("x",&x,hasBins); bindD("Q2",&Q2,hasBins); bindD("phi2",&phi2,hasPhi2); bindD("Delta_phi",&Delta_phi,hasDp);
        bindD("beam_pol",&beam_pol,hasPol);
        bindI("detector1",&detector1,hasTopo); bindI("detector2",&detector2,hasTopo);
    }
    double phi() const { return hasPhi2 ? phi2 : (hasDp ? Delta_phi : std::numeric_limits<double>::quiet_NaN()); }
};

static inline bool passesTopology_simple(int d1, int d2, const std::vector<std::string>& tops){
    for (const auto& t : tops){
        if (t=="(FD,FD)" && d1==1 && d2==1) return true;
        if (t=="(CD,FD)" && d1==2 && d2==1) return true;
        if (t=="(CD,FT)" && d1==2 && d2==0) return true;
    }
    return false;
}
static inline bool applyKinematicCuts_simple(double t1,double oa,double pT){ return !(oa<=5.0 || (-t1)>1.0 || pT>0.20); }

static std::vector<PolStats> compute_bin_polarization(
    TTree* t, const std::vector<std::pair<double,double>>& xB_bins,
    const std::vector<std::pair<double,double>>& Q2_bins,
    const std::vector<std::pair<double,double>>& t_bins,
    const std::vector<std::string>& tops)
{
    std::vector<PolStats> Pbin(xB_bins.size()*Q2_bins.size()*t_bins.size()*N_PHI_BINS);
    auto idx=[&](int ix,int iQ,int it,int ip){ return (((ix* (int)Q2_bins.size() + iQ)*(int)t_bins.size()+it)*N_PHI_BINS + ip); };

    BranchPol b; b.bind(t);
    if (!(b.hasH && b.hasCuts && b.hasBins && (b.hasPhi2||b.hasDp) && b.hasPol && b.hasTopo)) return Pbin;

    const Long64_t nent = t->GetEntries();
    for (Long64_t i=0;i<nent;++i){
        t->GetEntry(i);
        if (!applyKinematicCuts_simple(b.t1,b.open_angle_ep2,b.pTmiss)) continue;
        if (b.helicity!=+1 && b.helicity!=-1) continue;
        if (!passesTopology_simple(b.detector1,b.detector2,tops)) continue;

        double xB=b.x, Q2=b.Q2, tt=fabs(b.t1), phi=b.phi();
        if (!std::isfinite(xB)||!std::isfinite(Q2)||!std::isfinite(tt)||!std::isfinite(phi)) continue;

        auto findBin=[&](double v,const std::vector<std::pair<double,double>>& ranges)->int{
            for (int k=0;k<(int)ranges.size();++k) if (v>=ranges[k].first && v<ranges[k].second) return k;
            return -1;
        };
        int ix=findBin(xB,xB_bins), iQ=findBin(Q2,Q2_bins), it=findBin(tt,t_bins);
        if (ix<0||iQ<0||it<0) continue;

        double width = TWO_PI/double(N_PHI_BINS);
        double w = std::fmod(phi,TWO_PI); if (w<0) w+=TWO_PI; if (w>=TWO_PI) w = std::nextafter(TWO_PI,0.0);
        int ip = std::min(std::max(int(std::floor(w/width)),0), N_PHI_BINS-1);

        PolStats& ps = Pbin[idx(ix,iQ,it,ip)];
        ps.P += (b.beam_pol - ps.P)/(ps.n+1); // online mean
        ps.n++;
    }

    for (auto& ps : Pbin) ps.P_se = 0.0;
    return Pbin;
}

// ------------ math helpers for stabilizer penalties ------------
static inline double Dmin_over_domain(double B, double C) {
    // Let x = cosφ ∈ [-1,1]; D(φ) = 1 + B cosφ + C cos2φ
    // cos2φ = 2x^2 - 1 ⇒ D = 2C x^2 + B x + (1 - C)
    auto quad = [&](double x){ return 2.0*C*x*x + B*x + (1.0 - C); };
    double m = std::min(quad(-1.0), quad(+1.0)); // endpoints: 1 + C ∓ B
    if (C > 0.0) {
        const double x0 = -B/(4.0*C);
        if (x0 >= -1.0 && x0 <= 1.0) {
            const double fv = 1.0 - C - (B*B)/(8.0*C); // value at extremum inside
            m = std::min(m, fv);
        }
    }
    return m;
}
static inline double overBox(double A){ double v = std::fabs(A) - A_MAX_AMP; return (v>0.0)? v*v : 0.0; }

// ------------ BSA computation core ------------
struct BSApt { double phi=0.0; double bsa=0.0; double err=0.0; bool valid=false; };
struct FitRes { double A=0, Aerr=0, B1=0, B1err=0, B2=0, B2err=0, C=0, Cerr=0, chi2=0; int ndf=0; int status=0; };
struct CellResult {
    std::vector<BSApt> points; FitRes fit;
    double P_used=1.0; bool P_per_bin=false; double P_period_avg=1.0;
};

// Custom χ² + penalties minimization (Minuit2)
static FitRes fit_cell(const std::vector<BSApt>& pts){
    std::vector<double> x, y, ey;
    x.reserve(pts.size()); y.reserve(pts.size()); ey.reserve(pts.size());
    for (const auto& p : pts) {
        if (!p.valid) continue;
        x.push_back(p.phi);                 // radians
        y.push_back(p.bsa);
        ey.push_back(std::max(p.err, 1e-6));
    }
    FitRes fr; 
    const int n = (int)x.size();
    if (n < 4) { fr.status = 1; fr.ndf = 0; return fr; }

    // χ² + stabilizers
    auto chi2pen = [&](const double *par){
        const double A  = par[0];
        const double B1 = par[1];
        const double B2 = par[2];
        const double C  = par[3];

        double chi2 = 0.0;
        for (int i=0;i<n;++i){
            const double phi = x[i];
            double denom = 1.0 + B1*std::cos(phi) + B2*std::cos(2.0*phi);
            if (denom < EPS_DEN_EVAL) denom = EPS_DEN_EVAL; // evaluation clamp
            const double yhat = C + (A*std::sin(phi))/denom;
            const double pull = (y[i] - yhat)/ey[i];
            chi2 += pull*pull;
        }

        // (1) Denominator soft barrier over full φ-domain
        const double Dmin = Dmin_over_domain(B1, B2);
        double pen_den = 0.0;
        if (Dmin < EPS_DEN_FLOOR) {
            const double deficit = EPS_DEN_FLOOR - Dmin;
            pen_den = LAMBDA_DEN * deficit * deficit;
        }

        // (2) Amplitude box on denominator coefficients (and optional on A,C)
        double pen_amp = 0.0;
        pen_amp += overBox(B1);
        pen_amp += overBox(B2);
        // Optionally also constrain |A|,|C| very softly:
        // pen_amp += 0.25*overBox(A); pen_amp += 0.25*overBox(C);

        return chi2 + pen_den + LAMBDA_AMP*pen_amp;
    };

    // Minimizer
    std::unique_ptr<ROOT::Math::Minimizer> min(
        ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"));
    min->SetMaxFunctionCalls(10000);
    min->SetMaxIterations(10000);
    min->SetTolerance(1e-6);
    min->SetPrintLevel(0);

    ROOT::Math::Functor fcn(chi2pen, 4);
    min->SetFunction(fcn);

    // Parameters: A, B1, B2, C
    // Seeds kept modest; bounds keep things physical and help stability.
    min->SetLimitedVariable(0, "A",  0.10, 0.02, -1.0,  1.0);
    min->SetLimitedVariable(1, "B1", 0.00, 0.02, -1.0,  1.0);
    min->SetLimitedVariable(2, "B2", 0.00, 0.02, -1.0,  1.0);
    min->SetLimitedVariable(3, "C",  0.00, 0.02, -0.5,  0.5);

    const bool ok = min->Minimize();
    fr.status = ok ? 0 : 1;

    const double *par = min->X();
    const double *err = min->Errors();

    fr.A  = par[0]; fr.Aerr  = err[0];
    fr.B1 = par[1]; fr.B1err = err[1];
    fr.B2 = par[2]; fr.B2err = err[2];
    fr.C  = par[3]; fr.Cerr  = err[3];

    fr.chi2 = chi2pen(par);
    fr.ndf  = std::max(0, n - 4);
    return fr;
}

// JSON writers (unchanged)
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
              "\"B2\": {\"val\": "<<cr.fit.B2<<", \"err\": "<<cr.fit.B2err<<"}, "
              "\"C\": {\"val\": "<<cr.fit.C<<", \"err\": "<<cr.fit.Cerr<<"}, "
              "\"chi2\": "<<cr.fit.chi2<<", \"ndf\": "<<cr.fit.ndf<<", \"status\": "<<cr.fit.status<<"},\n";
        ofs<<"      \"polarization\": {\"per_bin\": "<<(cr.P_per_bin?"true":"false")<<", \"P_used\": "<<cr.P_used<<"}\n";
        ofs<<"    }";
    }
    ofs<<"\n  }\n}\n";
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
                  "\"B2\":{\"val\":"<<cr.fit.B2<<",\"err\":"<<cr.fit.B2err<<"},"
                  "\"C\":{\"val\":"<<cr.fit.C<<",\"err\":"<<cr.fit.Cerr<<"},"
                  "\"chi2\":"<<cr.fit.chi2<<",\"ndf\":"<<cr.fit.ndf<<",\"status\":"<<cr.fit.status<<"},";
            ofs<<"\"polarization\":{\"per_bin\":"<<(cr.P_per_bin?"true":"false")<<",\"P_used\":"<<cr.P_used<<"}";
            ofs<<"}";
        }
        ofs<<"\n      }\n    }";
    }
    ofs<<"\n  }\n}\n";
}

// ------------ plotting (dynamic grid + clean title + red fit curve) ------------
static void drawDegreeTicks(double xmin, double ymin, double xmax, double labelSize){
    // Draw 0,90,180,270,360 degree ticks with labels on top of the frame axis
    TGaxis* ax = new TGaxis(xmin, ymin, xmax, ymin, 0.0, 360.0, 4, "");
    ax->SetLabelFont(42);
    ax->SetLabelSize(labelSize);
    ax->SetLabelOffset(0.012);
    ax->SetTitle("");
    ax->SetTickSize(0.02);
    ax->Draw();
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
    using std::filesystem::create_directories;
    std::error_code ec;
    create_directories(out_dir_plots, ec);

    static const auto PHI_DEG = phiCentersDeg();

    for (int ix = 0; ix < (int)xB_bins.size(); ++ix) {
        const auto xb = xB_bins[ix];

        // Only the Q² and |t| bins that exist in this xB slice
        std::vector<std::pair<double,double>> Q2_slice, t_slice;
        uniqueQT_for_xB(binning_scheme, xb, Q2_slice, t_slice);
        if (Q2_slice.empty() || t_slice.empty()) continue;

        const int nrows = (int)t_slice.size();
        const int ncols = (int)Q2_slice.size();

        // Separate title pad (smaller now)
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

        // Title (compact and safe)
        pTop->cd();
        TLatex head;
        head.SetNDC(); head.SetTextAlign(22);
        head.SetTextFont(42);
        head.SetTextSize(0.36); // smaller than before
        std::ostringstream tit;
        tit << Form("Beam-Spin Asymmetry  %s   x_{B} #in [%.2g, %.2g]",
                    period.c_str(), xb.first, xb.second);
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

                // frame
                TH1* frame = gPad->DrawFrame(0.0, -1.05, 360.0, 1.05);
                TAxis* ax = frame->GetXaxis();
                TAxis* ay = frame->GetYaxis();

                // We will overlay a TGaxis for 0,90,180,270,360, so hide the default labels
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

                // overlay the custom degree ticks
                drawDegreeTicks(gPad->GetUxmin(), gPad->GetUymin(), gPad->GetUxmax(), 0.048);

                auto itCell = cells.find(std::make_tuple(ix, iQ_global, it_global));
                if (itCell == cells.end()) continue;
                const auto& cr = itCell->second;

                // Data points
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

                // Fitted curve (thin red line)
                if (cr.fit.status == 0 || cr.fit.ndf > 0) {
                    // draw via dense sampling to avoid any eval issues
                    const int NS=721;
                    std::vector<double> xd(NS), yd(NS);
                    for (int i=0;i<NS;++i){
                        double deg = double(i)*0.5;     // 0..360 in 0.5°
                        double rad = deg * (TWO_PI/360.0);
                        double denom = 1.0 + cr.fit.B1*std::cos(rad) + cr.fit.B2*std::cos(2.0*rad);
                        if (denom < EPS_DEN_EVAL) denom = EPS_DEN_EVAL;
                        double val = cr.fit.C + (cr.fit.A*std::sin(rad))/denom;
                        xd[i] = deg; yd[i] = val;
                    }
                    TGraph* gfit = new TGraph(NS, xd.data(), yd.data());
                    gfit->SetLineColor(kRed);
                    gfit->SetLineWidth(2);
                    gfit->Draw("L SAME");
                }

                // Panel annotation (Q² and -t)
                TLatex lab;
                lab.SetNDC(); lab.SetTextSize(0.040); lab.SetTextAlign(11);
                lab.SetTextFont(42);
                lab.DrawLatex(0.15, 0.94,
                    Form("Q^{2} #in [%.2g, %.2g],   -t #in [%.2g, %.2g]",
                         Q2_slice[ccol].first, Q2_slice[ccol].second,
                         t_slice[r].first,    t_slice[r].second));

                // Legend: keep well inside the pad
                TLegend* leg = new TLegend(0.60, 0.68, 0.88, 0.92);
                leg->SetBorderSize(1);
                leg->SetLineColor(kBlack);
                leg->SetFillStyle(0);
                leg->SetTextFont(42);
                leg->SetTextSize(0.040);
                leg->AddEntry((TObject*)nullptr, Form("A = %.3f #pm %.3f",  cr.fit.A,  cr.fit.Aerr), "");
                leg->AddEntry((TObject*)nullptr, Form("B_{1} = %.3f #pm %.3f", cr.fit.B1, cr.fit.B1err), "");
                leg->AddEntry((TObject*)nullptr, Form("B_{2} = %.3f #pm %.3f", cr.fit.B2, cr.fit.B2err), "");
                leg->AddEntry((TObject*)nullptr, Form("C = %.3f #pm %.3f",  cr.fit.C,  cr.fit.Cerr), "");
                leg->Draw();
            }
        }

        std::ostringstream fout;
        fout << out_dir_plots << "/plot_bsa_" << period << "_xB_" << ix << ".png";
        c->SaveAs(fout.str().c_str());
        delete c;
    }
}

// ------------ helpers ------------
static inline bool inTenSix(const std::string& runTag) {
    return (runTag=="sp18_inb" || runTag=="sp18_out" ||
            runTag=="fa18_inb" || runTag=="fa18_inb_supp" || runTag=="fa18_out");
}

} // anon namespace

// =======================================================
// Public driver
// =======================================================
void compute_and_plot_bsa_helicity(
    const std::vector<std::string>& periods,
    const std::vector<std::string>& topologies,
    const std::vector<Binning>& binning_scheme,
    const std::map<std::string, TTree*>& dvcsDataTrees,
    const std::string& total_counts_json_path,
    const std::string& contamination_dir,
    const std::string& out_root_dir)
{
    namespace fs = std::filesystem;

    const auto xB_bins = uniqueRanges(binning_scheme, 'x');
    const auto Q2_bins = uniqueRanges(binning_scheme, 'Q');
    const auto t_bins  = uniqueRanges(binning_scheme, 't');

    GroupCounts groups;
    if (!load_total_counts(total_counts_json_path, groups)) {
        std::cerr<<"[bsa][ERROR] Failed to load total_counts json.\n"; return;
    }

    const auto PHI_RAD = phiCentersRad();

    const fs::path json_period_dir = fs::path(out_root_dir)/"jsons"/"BSA_fits";
    fs::create_directories(json_period_dir);

    std::map<std::string, std::map<std::tuple<int,int,int>, CellResult>> allPeriodCells;

    struct Acc { double Np=0, Nm=0, Varp=0, Varm=0; };
    std::map<std::tuple<int,int,int,int>, Acc> acc106; // (ix,iQ,it,ip)

    for (const auto& period : periods) {
        const std::string runTag = periodToRunTagKey(period);

        auto itG = groups.find(runTag);
        if (itG == groups.end()) {
            std::cerr<<"[bsa][WARN] total_counts has no group '"<<runTag<<"'. Skipping "<<period<<"\n";
            continue;
        }
        const auto& countsMap = itG->second;

        ContamMap contam;
        const fs::path contam_path = fs::path(contamination_dir)/("contamination_"+period+".json");
        if (!load_contam_for_period(contam_path.string(), contam)) {
            std::cerr<<"[bsa][WARN] No contamination file for "<<period<<". Assuming c=0.\n";
        }

        // polarization
        TTree* t = nullptr;
        auto itT = dvcsDataTrees.find(runTag);
        if (itT!=dvcsDataTrees.end()) t = itT->second;
        std::vector<PolStats> Pbin;
        PolStats Pavg{1.0,0.0,0};
        if (t){
            Pbin = compute_bin_polarization(t, xB_bins, Q2_bins, t_bins, topologies);
            double sumP=0; int nn=0;
            for (auto& ps : Pbin) if (ps.n>0) { sumP+=ps.P; nn++; }
            if (nn>0) { Pavg.P = sumP/double(nn); Pavg.n=nn; }
        }

        std::map<std::tuple<int,int,int>, CellResult> cells;

        // build cells
        for (int ix=0; ix<(int)xB_bins.size(); ++ix)
        for (int iQ=0; iQ<(int)Q2_bins.size(); ++iQ)
        for (int itb=0; itb<(int)t_bins.size(); ++itb) {
            CellResult result;
            result.points.resize(N_PHI_BINS);
            result.P_per_bin = true;
            result.P_used = (Pavg.P>0? Pavg.P : 1.0);

            for (int ip=0; ip<N_PHI_BINS; ++ip){
                const BinKey bk(ix,iQ,itb,ip);

                // raw counts
                auto itC = countsMap.find(bk);
                long long Np=0, Nm=0;
                if (itC!=countsMap.end()) { Np = itC->second.plus; Nm = itC->second.minus; }

                // contamination
                ContamBin cb;
                auto itCt = contam.find(bk);
                if (itCt!=contam.end()) cb = itCt->second;

                // corrected counts
                double Np_corr = Np*(1.0 - cb.c_plus);
                double Nm_corr = Nm*(1.0 - cb.c_minus);
                double VarNp   = (1.0 - cb.c_plus)*(1.0 - cb.c_plus)*Np + (double(Np)*double(Np))*cb.c_plus_err*cb.c_plus_err;
                double VarNm   = (1.0 - cb.c_minus)*(1.0 - cb.c_minus)*Nm + (double(Nm)*double(Nm))*cb.c_minus_err*cb.c_minus_err;

                // polarization to use for this bin
                double P_here = result.P_used;
                if (!Pbin.empty()){
                    size_t idx = (((ix*(size_t)Q2_bins.size()+iQ)*(size_t)t_bins.size()+itb)*N_PHI_BINS + ip);
                    if (idx<Pbin.size() && Pbin[idx].n>0 && Pbin[idx].P>0.1) { P_here = Pbin[idx].P; }
                    else result.P_per_bin = false;
                }
                if (P_here<=0.0) { P_here = (Pavg.P>0? Pavg.P : 1.0); result.P_per_bin=false; }
                result.P_used = P_here;

                // scale by 1/P
                double Np_pol = Np_corr / P_here;
                double Nm_pol = Nm_corr / P_here;
                double VarNp_pol = VarNp / (P_here*P_here);
                double VarNm_pol = VarNm / (P_here*P_here);

                // BSA = (D/S)
                double S = Np_pol + Nm_pol;
                double D = Np_pol - Nm_pol;

                BSApt p; p.phi = PHI_RAD[ip];
                if (S>0.0) {
                    p.bsa = D/S;
                    double dNp = 2*Nm_pol/(S*S);
                    double dNm = -2*Np_pol/(S*S);
                    double VarBSA = dNp*dNp*VarNp_pol + dNm*dNm*VarNm_pol;
                    p.err = std::sqrt(std::max(VarBSA, 0.0));
                    p.valid = std::isfinite(p.bsa) && std::isfinite(p.err);
                } else {
                    p.bsa = 0.0; p.err = 0.0; p.valid=false;
                }
                result.points[ip]=p;

                // accumulate for combined 10.6
                if (inTenSix(runTag)) {
                    auto key4 = std::make_tuple(ix,iQ,itb,ip);
                    Acc& A = acc106[key4];
                    A.Np   += Np_pol;
                    A.Nm   += Nm_pol;
                    A.Varp += VarNp_pol;
                    A.Varm += VarNm_pol;
                }
            }

            result.fit = fit_cell(result.points);
            cells[std::make_tuple(ix,iQ,itb)] = std::move(result);
        }

        // write per-period JSON + plots
        const fs::path outP = json_period_dir/("BSA_fits_"+period+".json");
        write_period_bsa_json(outP.string(), N_PHI_BINS, xB_bins, Q2_bins, t_bins, cells);

        const fs::path plots_dir = fs::path(out_root_dir)/"bsa_plots"/periodToRunTagKey(period);
        std::error_code ec; fs::create_directories(plots_dir, ec);
        plot_cells_for_period(period, binning_scheme, xB_bins, Q2_bins, t_bins, cells, plots_dir.string());

        allPeriodCells[period] = std::move(cells);
    }

    // ---- write “all periods” file ----
    write_all_periods_json((fs::path(out_root_dir)/"jsons"/"BSA_fits_all_periods.json").string(),
                           allPeriodCells, N_PHI_BINS, xB_bins, Q2_bins, t_bins);

    // ---- 10.6 GeV statistical combination ----
    std::map<std::tuple<int,int,int>, CellResult> combCells;
    for (int ix=0; ix<(int)xB_bins.size(); ++ix)
    for (int iQ=0; iQ<(int)Q2_bins.size(); ++iQ)
    for (int itb=0; itb<(int)t_bins.size(); ++itb) {
        CellResult cr; cr.points.resize(N_PHI_BINS);
        cr.P_used = 1.0; cr.P_per_bin = false;
        for (int ip=0; ip<N_PHI_BINS; ++ip){
            auto itA = acc106.find(std::make_tuple(ix,iQ,itb,ip));
            BSApt p; p.phi = phiCentersRad()[ip]; p.valid=false;
            if (itA!=acc106.end()){
                const auto& A = itA->second;
                double S = A.Np + A.Nm;
                double D = A.Np - A.Nm;
                if (S>0.0) {
                    p.bsa = D/S;
                    double dNp = 2*A.Nm/(S*S);
                    double dNm = -2*A.Np/(S*S);
                    double VarBSA = dNp*dNp*A.Varp + dNm*dNm*A.Varm;
                    p.err = std::sqrt(std::max(VarBSA, 0.0));
                    p.valid = std::isfinite(p.bsa)&&std::isfinite(p.err);
                }
            }
            cr.points[ip]=p;
        }
        cr.fit = fit_cell(cr.points);
        combCells[std::make_tuple(ix,iQ,itb)] = std::move(cr);
    }

    // write combined JSON (uses 10.6 with a dot)
    {
        std::ofstream ofs((fs::path(out_root_dir)/"jsons"/"BSA_fits_combined_10.6.json").string());
        if (ofs){
            ofs<<std::fixed<<std::setprecision(8);
            ofs<<"{\n";
            ofs<<"  \"combined\": true,\n";
            ofs<<"  \"periods_used\": [\"DVCS_Sp18_inb\",\"DVCS_Sp18_out\",\"DVCS_Fa18_inb_supp\",\"DVCS_Fa18_inb\",\"DVCS_Fa18_out\"],\n";
            ofs<<"  \"binning_meta\": {\"phi_bins\": "<<N_PHI_BINS<<", \"xB_bins\": "<<xB_bins.size()<<", \"Q2_bins\": "<<Q2_bins.size()<<", \"t_bins\": "<<t_bins.size()<<"},\n";
            ofs<<"  \"bins\": {\n";
            bool first=true;
            for (const auto& kv : combCells){
                if (!first) ofs<<",\n"; first=false;
                int ix,iQ,it; std::tie(ix,iQ,it)=kv.first;
                const auto& cr = kv.second;
                ofs<<"    \"("<<ix<<","<<iQ<<","<<it<<")\": {";
                ofs<<"\"phi\":["; for (size_t i=0;i<cr.points.size();++i){ if(i)ofs<<","; ofs<<cr.points[i].phi; } ofs<<"],";
                ofs<<"\"bsa\":["; for (size_t i=0;i<cr.points.size();++i){ if(i)ofs<<","; ofs<<cr.points[i].bsa; } ofs<<"],";
                ofs<<"\"bsa_err\":["; for (size_t i=0;i<cr.points.size();++i){ if(i)ofs<<","; ofs<<cr.points[i].err; } ofs<<"],";
                ofs<<"\"fit\":{"
                      "\"A\":{\"val\":"<<cr.fit.A<<",\"err\":"<<cr.fit.Aerr<<"},"
                      "\"B1\":{\"val\":"<<cr.fit.B1<<",\"err\":"<<cr.fit.B1err<<"},"
                      "\"B2\":{\"val\":"<<cr.fit.B2<<",\"err\":"<<cr.fit.B2err<<"},"
                      "\"C\":{\"val\":"<<cr.fit.C<<",\"err\":"<<cr.fit.Cerr<<"},"
                      "\"chi2\":"<<cr.fit.chi2<<",\"ndf\":"<<cr.fit.ndf<<",\"status\":"<<cr.fit.status<<"}";
                ofs<<"}";
            }
            ofs<<"\n  }\n}\n";
        }
    }

    // plots for combined 10.6 (directory name matches the dot style)
    const fs::path plots_comb106 = fs::path(out_root_dir)/"bsa_plots"/"10.6_combined";
    std::error_code ec; fs::create_directories(plots_comb106, ec);
    plot_cells_for_period("RGA_10.6_combined", binning_scheme, xB_bins, Q2_bins, t_bins, combCells, plots_comb106.string());
}