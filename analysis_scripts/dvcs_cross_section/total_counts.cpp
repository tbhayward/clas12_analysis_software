// total_counts.cpp — produce per-period AND combined (Sp18, Fa18, 10.6_GeV) totals
// - Helicity-resolved DVCS data counts by (xB, Q2, |t|, phi)
// - Applies the same simple kinematic/topology cuts and optional 3σ cuts
// - JSON now contains a top-level "groups" object; each group has its own "bins"
// - Plots (multi-panel vs phi) for every group under output/total_counts_plots/<group>/...

#include "total_counts.h"

#include <TTree.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TH1.h>
#include <TGraphErrors.h>
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

// -------- helpers (consistent with your other files) --------
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

// runTag "fa18_inb" -> "DVCS_Fa18_inb"
static std::string dvcsPeriodName(const std::string& runTag) {
    std::string cap = runTag;
    if (!cap.empty()) cap[0] = static_cast<char>(std::toupper(static_cast<unsigned char>(cap[0])));
    for (size_t i = 0; i + 1 < cap.size(); ++i) {
        if (cap[i] == '_' && i + 1 < cap.size()) {
            cap[i+1] = static_cast<char>(std::toupper(static_cast<unsigned char>(cap[i+1])));
        }
    }
    return std::string("DVCS_") + cap;
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

// Unique (min,max) bins from scheme
static std::vector<std::pair<double,double>> uniqueRanges(const std::vector<Binning>& scheme, char which) {
    std::set<std::pair<double,double>> s;
    for (const auto& b : scheme) {
        if (which == 'x') s.emplace(b.xBmin, b.xBmax);
        else if (which == 'Q') s.emplace(b.Q2min, b.Q2max);
        else if (which == 't') s.emplace(b.tmin,  b.tmax);
    }
    return std::vector<std::pair<double,double>>(s.begin(), s.end());
}
static int findBin(double v, const std::vector<std::pair<double,double>>& ranges) {
    for (int i = 0; i < static_cast<int>(ranges.size()); ++i)
        if (v >= ranges[i].first && v < ranges[i].second) return i;
    return -1;
}

// --------- 3σ cuts loader (same parser style as contamination) ---------
struct Stats { double mean=0.0; double std=0.0; };
using VarCutMap = std::map<std::string, Stats>;
using PeriodTopoCuts = std::map<std::string, VarCutMap>; // key: "DVCS_Fa18_inb_FD_FD"

static bool loadCombinedCuts(const std::string& path, PeriodTopoCuts& out) {
    std::ifstream ifs(path);
    if (!ifs) { std::cerr << "[total_counts][WARN] Cannot open cuts JSON: " << path << "\n"; return false; }
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
static bool passes3SigmaCuts(const VarCutMap& cuts, const std::map<std::string,double>& values) {
    for (const auto& kv : cuts) {
        auto it = values.find(kv.first);
        if (it == values.end()) continue;
        if (!within3Sigma(it->second, kv.second)) return false;
    }
    return true;
}

// -------- branches --------
struct BranchBinderDVCS {
    int detector1=0, detector2=0; bool has_d1=false, has_d2=false;
    int helicity=0; bool has_helicity=false;
    double t1=0, open_angle_ep2=0, pTmiss=0; bool has_t1=false, has_oa=false, has_pT=false;
    double x=0, Q2=0, phi2=0, Delta_phi=0; bool has_x=false, has_Q2=false, has_phi2=false, has_Dp=false;

    // extra exclusivity vars to apply 3σ if present
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
        if (has_Dp)  m["Delta_phi"]=Delta_phi;
        if (has_tgg) m["theta_gamma_gamma"]=theta_gamma_gamma;
        if (has_pT)  m["pTmiss"]=pTmiss;
        if (has_xF)  m["xF"]=xF;
        if (has_Em)  m["Emiss2"]=Emiss2;
        if (has_Mx2) m["Mx2"]=Mx2;
        if (has_Mx21) m["Mx2_1"]=Mx2_1;
        if (has_Mx22) m["Mx2_2"]=Mx2_2;
        return m;
    }
};

// ---- containers ----
using BinKey = std::tuple<int,int,int,int>; // (ix,iQ2,it,ip)
struct HelCounts { long long plus=0, minus=0; };

static inline std::string keyStr(const BinKey& k) {
    int a,b,c,d; std::tie(a,b,c,d)=k;
    std::ostringstream os; os<<"("<<a<<","<<b<<","<<c<<","<<d<<")";
    return os.str();
}

// ---- JSON writer (has top-level "groups") ----
static void write_total_counts_json(
    const std::string& out_path,
    const std::map<std::string, std::map<BinKey, HelCounts>>& groups,
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

    bool firstG=true;
    for (const auto& gkv : groups) {
        if (!firstG) ofs << ",\n";
        firstG=false;
        ofs << "    \"" << gkv.first << "\": {\n";
        ofs << "      \"bins\": {\n";
        bool firstB=true;
        for (const auto& kv : gkv.second) {
            if (!firstB) ofs << ",\n";
            firstB=false;
            const HelCounts& hc = kv.second;
            ofs << "        \"" << keyStr(kv.first) << "\": {"
                << "\"helicity\":{\"+1\":" << hc.plus << ",\"-1\":" << hc.minus << "},"
                << "\"total\":" << (hc.plus + hc.minus)
                << "}";
        }
        ofs << "\n      }\n    }";
    }
    ofs << "\n  }\n";
    ofs << "}\n";
    ofs.close();
    std::cout << "[total_counts] Wrote " << out_path << "\n";
}

// ---------------- plotting (multi-panel vs φ, rows=Q², cols=|t|) ----------------
static std::vector<double> phiCentersDeg(int nPhi=N_PHI_BINS) {
    std::vector<double> v(nPhi);
    const double step = 360.0 / static_cast<double>(nPhi);
    for (int i = 0; i < nPhi; ++i) v[i] = (i + 0.5) * step;
    return v;
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
static int findIndex(const std::pair<double,double>& range,
                     const std::vector<std::pair<double,double>>& ranges) {
    for (int i = 0; i < static_cast<int>(ranges.size()); ++i) {
        if (ranges[i] == range) return i;
    }
    return -1;
}

static void plot_group_counts(
    const std::string& group_name,
    const std::map<BinKey, HelCounts>& table,
    const std::vector<Binning>& binning_scheme,
    const std::vector<std::pair<double,double>>& xB_bins,
    const std::vector<std::pair<double,double>>& Q2_bins,
    const std::vector<std::pair<double,double>>& t_bins,
    const std::string& out_dir_plots
){
    using namespace std;
    namespace fs = std::filesystem;
    std::error_code ec;
    fs::create_directories(out_dir_plots, ec);

    static const std::vector<double> PHI = phiCentersDeg();
    std::vector<double> X(N_PHI_BINS), ex(N_PHI_BINS, 0.0);
    for (int i=0;i<N_PHI_BINS;++i) X[i] = PHI[i];

    for (int ix=0; ix<(int)xB_bins.size(); ++ix) {
        const auto xBr = xB_bins[ix];

        vector<pair<double,double>> Q2_slice, t_slice;
        uniqueQT_for_xB(binning_scheme, xBr, Q2_slice, t_slice);
        if (Q2_slice.empty() || t_slice.empty()) continue;

        const int nrows = (int)Q2_slice.size();
        const int ncols = (int)t_slice.size();

        const int W = 260*ncols + 120;
        const int H = 220*nrows + 120;

        ostringstream cname; cname<<"c_counts_"<<group_name<<"_xB"<<ix;
        TCanvas* c = new TCanvas(cname.str().c_str(), cname.str().c_str(), W, H);

        TPad* pTop  = new TPad("pTop","pTop", 0.0, 0.94, 1.0, 1.0);
        pTop->SetFillStyle(0); pTop->SetBorderSize(0); pTop->Draw();

        TPad* pGrid = new TPad("pGrid","pGrid", 0.0, 0.0, 1.0, 0.94);
        pGrid->SetFillStyle(0); pGrid->SetBorderSize(0); pGrid->Draw();
        pGrid->cd();
        pGrid->Divide(ncols, nrows, 0.0001, 0.0001);

        pTop->cd();
        TLatex head;
        head.SetNDC();
        head.SetTextSize(0.55);
        head.SetTextAlign(22);
        head.DrawLatex(0.5, 0.5, group_name.c_str());

        for (int r=0; r<nrows; ++r) {
            int iQ = findIndex(Q2_slice[r], Q2_bins); if (iQ<0) continue;
            for (int cc=0; cc<ncols; ++cc) {
                int itb = findIndex(t_slice[cc], t_bins); if (itb<0) continue;

                pGrid->cd(r*ncols + cc + 1);
                gPad->SetGrid(1,1);
                gPad->SetTopMargin(0.08);
                gPad->SetBottomMargin(0.14);
                gPad->SetLeftMargin(0.12);
                gPad->SetRightMargin(0.06);

                std::vector<double> Yp(N_PHI_BINS, 0.0), Ym(N_PHI_BINS, 0.0);
                std::vector<double> eYp(N_PHI_BINS, 0.0), eYm(N_PHI_BINS, 0.0);

                double ymax_found = 0.0;
                for (int ip=0; ip<N_PHI_BINS; ++ip) {
                    auto it = table.find(BinKey(ix, iQ, itb, ip));
                    if (it == table.end()) continue;
                    const HelCounts& hc = it->second;

                    Yp[ip]  = (double)hc.plus;
                    Ym[ip]  = (double)hc.minus;
                    eYp[ip] = std::sqrt(std::max(0.0, Yp[ip]));
                    eYm[ip] = std::sqrt(std::max(0.0, Ym[ip]));
                    ymax_found = std::max(ymax_found, std::max(Yp[ip] + eYp[ip], Ym[ip] + eYm[ip]));
                }

                const double ymin = 0.0;
                const double ymax = (ymax_found > 0.0 ? ymax_found*1.25 : 10.0);

                TH1 *frame = gPad->DrawFrame(0.0, ymin, 360.0, ymax);
                frame->GetXaxis()->SetTitle("#phi (deg)");
                frame->GetYaxis()->SetTitle("DVCS counts");
                frame->GetXaxis()->SetNdivisions(505);
                frame->GetXaxis()->CenterTitle();
                frame->GetYaxis()->CenterTitle();

                TGraphErrors* grP = new TGraphErrors(N_PHI_BINS, X.data(), Yp.data(), ex.data(), eYp.data());
                TGraphErrors* grM = new TGraphErrors(N_PHI_BINS, X.data(), Ym.data(), ex.data(), eYm.data());
                grP->SetMarkerStyle(24);
                grM->SetMarkerStyle(20);
                grP->SetMarkerSize(0.9);
                grM->SetMarkerSize(0.9);
                grP->SetLineWidth(2);
                grM->SetLineWidth(2);
                grP->Draw("P SAME");
                grM->Draw("P SAME");

                TLegend* leg = new TLegend(0.60, 0.73, 0.92, 0.92);
                leg->SetBorderSize(1);
                leg->SetLineColor(kBlack);
                leg->SetFillStyle(0);
                leg->SetTextSize(0.035);
                leg->AddEntry(grP, "helicity +1", "p");
                leg->AddEntry(grM, "helicity -1", "p");
                leg->Draw();
            }
        }

        ostringstream fout;
        fout<<out_dir_plots<<"/plot_counts_"<<group_name<<"_xB_"<<ix<<".png";
        c->SaveAs(fout.str().c_str());
        delete c;
    }
}

// ---- small helpers to build combined groups ----
static inline bool isSpring18(const std::string& period) {
    std::string s = toLower(period);
    return (s.find("sp18") != std::string::npos || s.find("spring18") != std::string::npos);
}
static inline bool isFall18(const std::string& period) {
    std::string s = toLower(period);
    return (s.find("fa18") != std::string::npos || s.find("fall18") != std::string::npos);
}
static void add_into(std::map<BinKey, HelCounts>& dst, const std::map<BinKey, HelCounts>& src) {
    for (const auto& kv : src) {
        auto& d = dst[kv.first];
        d.plus  += kv.second.plus;
        d.minus += kv.second.minus;
    }
}

// ---- main ----
// NOTE: signature changed to include out_root_dir for plots in a sibling folder to jsons.
void compute_total_counts(
    const std::vector<std::string>& periods,
    const std::vector<std::string>& topologies,
    const std::vector<Binning>& binning_scheme,
    const std::map<std::string, TTree*>& dvcsDataTrees,
    const std::string& combined_cuts_json,
    const std::string& out_json_path,
    const std::string& out_root_dir  // e.g. "output"
) {
    namespace fs = std::filesystem;

    // Build global bin edges
    const auto xB_bins = uniqueRanges(binning_scheme, 'x');
    const auto Q2_bins = uniqueRanges(binning_scheme, 'Q');
    const auto t_bins  = uniqueRanges(binning_scheme, 't');
    if (xB_bins.empty() || Q2_bins.empty() || t_bins.empty()) {
        std::cerr << "[total_counts][ERROR] Binning ranges are empty.\n";
        return;
    }

    // Load combined 3σ cuts (data)
    PeriodTopoCuts cuts;
    if (!loadCombinedCuts(combined_cuts_json, cuts)) {
        std::cerr << "[total_counts][WARN] No combined cuts loaded; proceeding without 3σ cuts.\n";
    }

    auto dvcsCutsKey = [&](const std::string& runTag, const std::string& topoKey)->std::string {
        std::string cap = runTag;
        if (!cap.empty()) cap[0] = std::toupper(static_cast<unsigned char>(cap[0]));
        for (size_t i=0;i+1<cap.size();++i) if (cap[i]=='_' && i+1<cap.size()) cap[i+1]=std::toupper(static_cast<unsigned char>(cap[i+1]));
        return std::string("DVCS_") + cap + "_" + topoKey;
    };

    // group name -> (bin -> counts)
    std::map<std::string, std::map<BinKey, HelCounts>> outCounts;

    // Per-period filling
    for (const auto& period : periods) {
        const std::string runTag = periodToRunTagKey(period);  // e.g. "fa18_inb"
        const std::string dvcsKeyInMap = runTag;               // tree key in dvcsDataTrees

        auto it = dvcsDataTrees.find(dvcsKeyInMap);
        if (it == dvcsDataTrees.end() || !it->second) {
            std::cerr << "[total_counts][WARN] No DVCS data tree for '" << period
                      << "' (key '" << dvcsKeyInMap << "'). Skipping.\n";
            continue;
        }
        TTree* t = it->second;

        BranchBinderDVCS b; b.bind(t);
        if (!b.readyCuts() || !b.readyBins()) {
            std::cerr << "[total_counts][WARN] DVCS tree missing branches for '" << period << "'. Skipping.\n";
            continue;
        }

        const Long64_t nent = t->GetEntries();
        for (Long64_t i=0;i<nent;++i) {
            t->GetEntry(i);
            if (!applyKinematicCuts_simple(b.t1,b.open_angle_ep2,b.pTmiss)) continue;
            if (b.helicity!=+1 && b.helicity!=-1) continue;

            // choose topology
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
                // collect values for cut variables
                std::map<std::string,double> vals;
                vals["Delta_phi"] = b.Delta_phi;
                vals["theta_gamma_gamma"] = b.theta_gamma_gamma;
                vals["pTmiss"] = b.pTmiss;
                vals["xF"] = b.xF;
                vals["Emiss2"] = b.Emiss2;
                vals["Mx2"] = b.Mx2;
                vals["Mx2_1"] = b.Mx2_1;
                vals["Mx2_2"] = b.Mx2_2;
                if (!passes3SigmaCuts(topoCuts, vals)) continue;
            }

            // bin
            double xB=b.x, Q2=b.Q2, tt=std::fabs(b.t1), phi=(b.phi2 || b.Delta_phi ? (b.phi2 ? b.phi2 : b.Delta_phi) : std::numeric_limits<double>::quiet_NaN());
            if (!std::isfinite(xB)||!std::isfinite(Q2)||!std::isfinite(tt)||!std::isfinite(phi)) continue;
            int ix=findBin(xB,xB_bins), iQ=findBin(Q2,Q2_bins), itb=findBin(tt,t_bins), ip=phiToBin(phi);
            if (ix<0||iQ<0||itb<0||ip<0||ip>=N_PHI_BINS) continue;
            BinKey key(ix,iQ,itb,ip);

            // (1) plain runTag (what downstream expects)
            HelCounts& hc_plain = outCounts[runTag][key];
            if (b.helicity==+1) hc_plain.plus++; else hc_plain.minus++;

            // (2) ALSO under DVCS_* name
            const std::string dvcsName = dvcsPeriodName(runTag);
            HelCounts& hc_dvcs = outCounts[dvcsName][key];
            if (b.helicity==+1) hc_dvcs.plus++; else hc_dvcs.minus++;

            // We’ll create combined groups after the per-period loops
        }
    }

    // ----- build the combined groups -----
    // Spring 2018 (sum all periods whose name contains Sp18)
    std::map<BinKey, HelCounts> spring_sum, fall_sum, all_sum;

    for (const auto& period : periods) {
        const std::string runTag = periodToRunTagKey(period);
        auto it = outCounts.find(runTag);
        if (it == outCounts.end()) continue;
        if (isSpring18(period)) add_into(spring_sum, it->second);
        if (isFall18(period))   add_into(fall_sum,   it->second);
        add_into(all_sum, it->second);
    }

    if (!spring_sum.empty()) outCounts["Spring2018"] = std::move(spring_sum);
    if (!fall_sum.empty())   outCounts["Fall2018"]   = std::move(fall_sum);
    if (!all_sum.empty())    outCounts["10.6_GeV"]   = std::move(all_sum); // exact name requested

    // ----- write JSON (with top-level "groups") -----
    write_total_counts_json(out_json_path, outCounts, N_PHI_BINS,
                            xB_bins, Q2_bins, t_bins);

    // ----- plots for every group -----
    const fs::path plots_root = fs::path(out_root_dir) / "total_counts_plots";
    for (const auto& gkv : outCounts) {
        const std::string& group = gkv.first;
        const fs::path out_dir = plots_root / group;
        plot_group_counts(group, gkv.second, binning_scheme, xB_bins, Q2_bins, t_bins, out_dir.string());
    }
}