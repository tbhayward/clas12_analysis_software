// exclusivity_cuts.cpp
// Same plots; universal cuts use global_cuts.h so every stage can share them.

#include "exclusivity_cuts.h"
#include "periods.h"
#include "global_cuts.h"

// ROOT
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TF1.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TPaveText.h>
#include <TROOT.h>
#include <TH1.h>
#include <TString.h>
#include <TDataType.h>  // kInt_t, kDouble_t and EDataType
#include <TSystem.h>
#include <TLine.h>

// C++ stdlib
#include <algorithm>
#include <atomic>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <mutex>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <set>
#include <thread>
#include <chrono>
#include <utility>
#include <vector>

static constexpr double PI = 3.14159265358979323846;

static inline double delta_phi_rad_from_two_phi(double phi_a, double phi_b) {
    double d = std::fmod(phi_a - phi_b, 2.0 * PI);

    if (d <= -PI) {
        d += 2.0 * PI;
    }

    if (d > PI) {
        d -= 2.0 * PI;
    }

    return std::fabs(d);
}

// -------------------- helpers: strings and keys --------------------

static std::string channelToStr(Channel ch) {
    return (ch == Channel::DVCS) ? "dvcs" : "eppi0";
}

static std::string channelPretty(Channel ch) {
    return (ch == Channel::DVCS) ? "DVCS" : "eppi0";
}

static std::string topoToKey(Topology t) {
    switch (t) {
        case Topology::FD_FD: return "FD_FD";
        case Topology::CD_FD: return "CD_FD";
        case Topology::CD_FT: return "CD_FT";
    }
    return "FD_FD";
}

static std::string periodCode(Channel ch, const std::string& runTag) {
    std::string prefix = (ch == Channel::DVCS) ? "DVCS_" : "eppi0_";
    std::string nice = runTag;
    if (!nice.empty()) nice[0] = std::toupper(nice[0]);
    for (size_t i = 0; i + 1 < nice.size(); ++i) {
        if (nice[i] == '_' && i + 1 < nice.size()) nice[i + 1] = std::toupper(nice[i + 1]);
    }
    return prefix + nice;
}

// -------------------- histogram configs (ordered) --------------------

using HistList = std::vector<std::pair<std::string, HistCfg>>;

static HistList getHistConfigs(Channel ch) {
    HistList cfg;
    cfg.push_back({"Delta_phi",            {100, 2.84159, 3.44159}});
    if (ch == Channel::DVCS) cfg.push_back({"theta_gamma_gamma", {100, 0.0, 2.0}});
    else                     cfg.push_back({"theta_pi0_pi0",     {100, 0.0, 2.0}});
    cfg.push_back({"pTmiss",              {100, 0.0, 0.3}});
    cfg.push_back({"xF",                  {100, -0.4, 0.2}});
    cfg.push_back({"Emiss2",              {100, -1.0, 2.0}});
    cfg.push_back({"Mx2",                 {100, -0.03, 0.03}});
    cfg.push_back({"Mx2_1",               {100, -1.5, 1.5}});
    cfg.push_back({"Mx2_2",               {100, 0.0, 3.0}});
    return cfg;
}

// -------------------- stage plan --------------------

static std::vector<StagePlan> buildStages(Channel ch) {
    std::vector<StagePlan> stages;
    stages.push_back(StagePlan{{"Mx2", "Mx2_1"}});       // stage 0
    stages.push_back(StagePlan{{"Emiss2", "Mx2_2"}});    // stage 1
    if (ch == Channel::DVCS) stages.push_back(StagePlan{{"Delta_phi", "pTmiss", "xF", "theta_gamma_gamma"}});
    else                     stages.push_back(StagePlan{{"Delta_phi", "pTmiss", "xF", "theta_pi0_pi0"}});
    return stages; // final pass is len(stages)+1 in driver
}

// -------------------- topology --------------------

static bool passesTopology(int detector1, int detector2, Topology topo) {
    if (topo == Topology::FD_FD) return (detector1 == 1 && detector2 == 1);
    if (topo == Topology::CD_FD) return (detector1 == 2 && detector2 == 1);
    if (topo == Topology::CD_FT) return (detector1 == 2 && detector2 == 0);
    return false;
}

static bool isUpperTailQuantileVar(const std::string& var) {
    return (var == "pTmiss" || var == "theta_gamma_gamma" || var == "theta_pi0_pi0");
}

static bool isUpperTailCut(const Stats& s) {
    return (s.mode == "upper_quantile");
}

static bool withinCutWindow(double val, const Stats& s) {
    if (!std::isfinite(val)) return false;

    if (isUpperTailCut(s)) {
        if (!std::isfinite(s.cut_high)) return true;
        return val <= s.cut_high;
    }

    double lo = s.cut_low;
    double hi = s.cut_high;

    if (!(std::isfinite(lo) && std::isfinite(hi)) || hi <= lo) {
        if (!std::isfinite(s.mean) || !std::isfinite(s.std) || s.std <= 0.0) return true;
        lo = s.mean - 3.0 * s.std;
        hi = s.mean + 3.0 * s.std;
    }

    return (val >= lo && val <= hi);
}

static bool passes3SigmaCuts(const std::map<std::string, Stats>& cuts,
                             const std::map<std::string, double>& values) {
    for (const auto& kv : cuts) {
        auto it = values.find(kv.first);
        if (it == values.end()) continue;
        if (!withinCutWindow(it->second, kv.second)) return false;
    }
    return true;
}

// -------------------- branch binder --------------------

static std::mutex g_root_bind_mutex;
static std::mutex g_fit_mutex;

struct BranchBinder {
    int runnum = 0;       bool has_runnum = false;
    int detector1 = 0;    bool has_detector1 = false;
    int detector2 = 0;    bool has_detector2 = false;

    double t1 = 0.0;                bool has_t1 = false;
    double open_angle_ep2 = 0.0;    bool has_open_angle_ep2 = false;
    double Emiss2 = 0.0;            bool has_Emiss2 = false;
    double Mx2 = 0.0;               bool has_Mx2 = false;
    double Mx2_1 = 0.0;             bool has_Mx2_1 = false;
    double Mx2_2 = 0.0;             bool has_Mx2_2 = false;
    double pTmiss = 0.0;            bool has_pTmiss = false;
    double xF = 0.0;                bool has_xF = false;
    double Delta_phi = 0.0;         bool has_Delta_phi = false;
    double theta_gamma_gamma = 0.0; bool has_theta_gamma_gamma = false;
    double theta_pi0_pi0 = 0.0;     bool has_theta_pi0_pi0 = false;

    // NEW (only used if cfg.enable_dvcsgen_ycol_cut is true): scattered e- and photon kinematics
    double e_p = 0.0;       bool has_e_p = false;
    double e_theta = 0.0;   bool has_e_theta = false;
    double e_phi = 0.0;     bool has_e_phi = false;

    double p1_theta = 0.0;  bool has_p1_theta = false;
    double p1_phi = 0.0;    bool has_p1_phi = false;

    double p2_p = 0.0;      bool has_p2_p = false;
    double p2_theta = 0.0;  bool has_p2_theta = false;
    double p2_phi = 0.0;    bool has_p2_phi = false;

    void bind(TTree* t, Channel ch) {
        if (!t) return;

        // ROOT branch-status/address operations touch shared ROOT internals.
        // Serialize binding while still allowing independent event loops to run
        // concurrently after binding is complete.
        std::lock_guard<std::mutex> lock(g_root_bind_mutex);

        // Disable everything, then explicitly enable only what we use.
        t->SetBranchStatus("*", 0);

        auto ena = [&](const char* n){ if (t->GetBranch(n)) t->SetBranchStatus(n, 1); };

        ena("runnum");
        ena("detector1");
        ena("detector2");
        ena("t1");
        ena("open_angle_ep2");
        ena("Emiss2");
        ena("Mx2");
        ena("Mx2_1");
        ena("Mx2_2");
        ena("pTmiss");
        ena("xF");
        ena("Delta_phi");
        if (ch == Channel::DVCS) ena("theta_gamma_gamma");
        else                     ena("theta_pi0_pi0");

        // NEW: required only when dvcsgen ycol cut is enabled.
        ena("e_p");
        ena("e_theta");
        ena("e_phi");
        ena("p1_theta");
        ena("p1_phi");
        ena("p2_p");
        ena("p2_theta");
        ena("p2_phi");

        // Bind with explicit primitive types (6-arg overload).
        auto bI = [&](const char* n, int* a, bool& f){
            if (t->GetBranch(n)) { t->SetBranchAddress(n, a, nullptr, nullptr, kInt_t, false); f = true; }
        };
        auto bD = [&](const char* n, double* a, bool& f){
            if (t->GetBranch(n)) { t->SetBranchAddress(n, a, nullptr, nullptr, kDouble_t, false); f = true; }
        };

        bI("runnum",    &runnum,    has_runnum);
        bI("detector1", &detector1, has_detector1);
        bI("detector2", &detector2, has_detector2);

        bD("t1", &t1, has_t1);
        bD("open_angle_ep2", &open_angle_ep2, has_open_angle_ep2);

        bD("Emiss2", &Emiss2, has_Emiss2);
        bD("Mx2", &Mx2, has_Mx2);
        bD("Mx2_1", &Mx2_1, has_Mx2_1);
        bD("Mx2_2", &Mx2_2, has_Mx2_2);
        bD("pTmiss", &pTmiss, has_pTmiss);
        bD("xF", &xF, has_xF);
        bD("Delta_phi", &Delta_phi, has_Delta_phi);
        if (ch == Channel::DVCS) bD("theta_gamma_gamma", &theta_gamma_gamma, has_theta_gamma_gamma);
        else                     bD("theta_pi0_pi0",     &theta_pi0_pi0,     has_theta_pi0_pi0);

        // NEW (may or may not exist, depending on tree schema)
        bD("e_p",     &e_p,     has_e_p);
        bD("e_theta", &e_theta, has_e_theta);
        bD("e_phi",   &e_phi,   has_e_phi);
        bD("p1_theta", &p1_theta, has_p1_theta);
        bD("p1_phi",  &p1_phi,  has_p1_phi);

        bD("p2_p",     &p2_p,     has_p2_p);
        bD("p2_theta", &p2_theta, has_p2_theta);
        bD("p2_phi",   &p2_phi,   has_p2_phi);
    }

    double delta_phi_value(bool& has_val) const {
        if (has_Delta_phi) {
            has_val = true;
            return Delta_phi;
        }

        if (has_p1_phi && has_p2_phi) {
            has_val = true;
            return delta_phi_rad_from_two_phi(p1_phi, p2_phi);
        }

        has_val = false;
        return 0.0;
    }

    std::map<std::string, double> valuesMap(Channel ch) const {
        std::map<std::string, double> m;
        bool has_delta_phi = false;
        const double delta_phi = delta_phi_value(has_delta_phi);
        if (has_delta_phi) m["Delta_phi"] = delta_phi;
        if (ch == Channel::DVCS) { if (has_theta_gamma_gamma) m["theta_gamma_gamma"] = theta_gamma_gamma; }
        else { if (has_theta_pi0_pi0) m["theta_pi0_pi0"] = theta_pi0_pi0; }
        if (has_pTmiss) m["pTmiss"] = pTmiss;
        if (has_xF) m["xF"] = xF;
        if (has_Emiss2) m["Emiss2"] = Emiss2;
        if (has_Mx2) m["Mx2"] = Mx2;
        if (has_Mx2_1) m["Mx2_1"] = Mx2_1;
        if (has_Mx2_2) m["Mx2_2"] = Mx2_2;
        return m;
    }
};

// -------------------- fits and stats --------------------

static std::pair<double,double> fitGaussianLeftSide(TH1D* h) {
    if (!h || h->GetEntries() == 0) return {0.0, 0.0};

    // ROOT fitting creates named objects and touches global fitter state.
    // Serialize this small section to keep the surrounding period-level
    // parallelism safe and deterministic.
    std::lock_guard<std::mutex> lock(g_fit_mutex);

    double xmin = h->GetXaxis()->GetXmin();
    int peakBin = h->GetMaximumBin();
    double xpeak = h->GetXaxis()->GetBinCenter(peakBin);
    double xr = std::max(xmin + 1e-6, 1.35 * xpeak);

    const std::string fname = std::string("gausLeft_") + h->GetName();
    TF1 f(fname.c_str(), "gaus(0)", xmin, xr);
    f.SetParameters(h->GetMaximum(), h->GetMean(), std::max(1e-3, h->GetRMS() / 2.0));
    f.SetParLimits(2, 1e-6, 5.0);
    int status = h->Fit(&f, "RQ0");
    if (status != 0) return {h->GetMean(), h->GetRMS()};
    double mu = f.GetParameter(1);
    double sg = std::abs(f.GetParameter(2));
    if (sg <= 1e-9) sg = h->GetRMS();
    return {mu, sg};
}

static Stats meanStd(TH1D* h) {
    if (!h || h->GetEntries() == 0) return {0.0, 0.0};
    return {h->GetMean(), h->GetStdDev()};
}

static void normalizeHist(TH1D* h) {
    if (!h) return;
    double integral = h->Integral();
    if (integral > 0.0) h->Scale(1.0 / integral);
}

// -------------------- pretty x-axis labels --------------------

static std::string formatLabelName(const std::string& var, Channel ch) {
    if (var == "Delta_phi")          return "#Delta#phi (rad)";
    if (var == "theta_gamma_gamma")  return "#theta_{#gamma#gamma} (rad)";
    if (var == "theta_pi0_pi0")      return "#theta_{#pi_{0}#pi_{0}} (rad)";
    if (var == "pTmiss")             return "p_{T}^{miss} (GeV)";
    if (var == "xF")                 return "x_{F}";
    if (var == "Emiss2")             return "E_{miss}^{2} (GeV^{2})";
    if (var == "Mx2")                return "M_{x}^{2} (GeV^{2})";
    if (var == "Mx2_1")              return "M_{x1}^{2} (GeV^{2})";
    if (var == "Mx2_2")              return "M_{x2}^{2} (GeV^{2})";
    return var;
}

// -------------------- stage fill --------------------

struct FilledHists {
    std::map<std::string, TH1D*> data;
    std::map<std::string, TH1D*> mc;
};

using FilledHistsByTopology = std::map<Topology, FilledHists>;

static Topology topologyFromDetectors(int detector1, int detector2, bool& valid) {
    valid = true;
    if (detector1 == 1 && detector2 == 1) return Topology::FD_FD;
    if (detector1 == 2 && detector2 == 1) return Topology::CD_FD;
    if (detector1 == 2 && detector2 == 0) return Topology::CD_FT;
    valid = false;
    return Topology::FD_FD;
}

static GlobalCutConfig globalConfigForStageTopology(Topology topo) {
    GlobalCutConfig gcfg = default_global_cuts();
    gcfg.enable_topology_filter = true;

    if (topo == Topology::FD_FD) {
        gcfg.required_detector1 = 1;
        gcfg.required_detector2 = 1;
    } else if (topo == Topology::CD_FD) {
        gcfg.required_detector1 = 2;
        gcfg.required_detector2 = 1;
    } else {
        gcfg.required_detector1 = 2;
        gcfg.required_detector2 = 0;
    }

    return gcfg;
}

static bool passesGlobalForStage(const BranchBinder& b,
                                 const std::string& period_label,
                                 const GlobalCutConfig& gcfg) {
    if (!(b.has_t1 && b.has_open_angle_ep2)) return false;
    if (gcfg.enable_pTmiss_cut && !b.has_pTmiss) return false;
    if (!(b.has_detector1 && b.has_detector2)) return false;

    if (global_cuts_require_sector_phi(gcfg)) {
        if (!(b.has_e_phi && b.has_p1_phi && b.has_p2_phi)) {
            throw std::runtime_error("[exclusivity_cuts] FATAL: sector selection requires e_phi, p1_phi, p2_phi.");
        }
    }

    if (global_cuts_require_auxiliary_kinematics(gcfg)) {
        if (!(b.has_e_theta && b.has_e_phi &&
              b.has_p1_theta && b.has_p1_phi &&
              b.has_p2_p && b.has_p2_theta && b.has_p2_phi)) {
            throw std::runtime_error("[exclusivity_cuts] FATAL: auxiliary fiducial cuts require e_theta, e_phi, p1_theta, p1_phi, p2_p, p2_theta, p2_phi.");
        }

        return passes_global_cuts(b.t1, b.open_angle_ep2, b.pTmiss,
                                  b.detector1, b.detector2,
                                  period_label,
                                  b.e_p, b.e_theta, b.e_phi,
                                  b.p1_theta, b.p1_phi,
                                  b.p2_p, b.p2_theta, b.p2_phi,
                                  gcfg);
    }

    if (gcfg.enable_dvcsgen_ycol_cut) {
        if (!(b.has_e_p && b.has_e_theta && b.has_e_phi &&
              b.has_p2_p && b.has_p2_theta && b.has_p2_phi)) {
            throw std::runtime_error("[exclusivity_cuts] FATAL: dvcsgen ycol cut requires e_p, e_theta, e_phi, p2_p, p2_theta, p2_phi.");
        }

        if (global_cuts_require_sector_phi(gcfg)) {
            return passes_global_cuts(b.t1, b.open_angle_ep2, b.pTmiss,
                                      b.detector1, b.detector2,
                                      period_label,
                                      b.e_p, b.e_theta, b.e_phi, b.p1_theta, b.p1_phi,
                                      b.p2_p, b.p2_theta, b.p2_phi,
                                      gcfg);
        }

        return passes_global_cuts(b.t1, b.open_angle_ep2, b.pTmiss,
                                  b.detector1, b.detector2,
                                  period_label,
                                  b.e_p, b.e_theta, b.e_phi,
                                  b.p2_p, b.p2_theta, b.p2_phi,
                                  gcfg);
    }

    if (global_cuts_require_sector_phi(gcfg)) {
        return passes_global_cuts(b.t1, b.open_angle_ep2, b.pTmiss,
                                  b.detector1, b.detector2,
                                  b.e_phi, b.p1_phi, b.p2_phi,
                                  gcfg);
    }

    return passes_global_cuts(b.t1, b.open_angle_ep2, b.pTmiss,
                              b.detector1, b.detector2,
                              gcfg);
}

static FilledHistsByTopology fillStageHistsAllTopologies(
    TTree* dataTree, TTree* mcTree, Channel ch,
    const std::string& period_label,
    const std::map<Topology, CutDict>& cumulative,
    const HistList& cfg, int stage_index)
{
    static const Topology topologies[] = {
        Topology::FD_FD, Topology::CD_FD, Topology::CD_FT
    };

    FilledHistsByTopology out;
    std::map<Topology, GlobalCutConfig> globalConfigs;

    for (Topology topo : topologies) {
        FilledHists& hists = out[topo];
        globalConfigs.emplace(topo, globalConfigForStageTopology(topo));

        const std::string histTag = channelToStr(ch) + "_" + period_label + "_" +
                                    topoToKey(topo) + "_stage" + std::to_string(stage_index);

        for (const auto& kv : cfg) {
            const std::string& var = kv.first;
            const HistCfg& hc = kv.second;

            auto* dh = new TH1D(("data_" + histTag + "_" + var).c_str(), "",
                                hc.nbins, hc.xlow, hc.xhigh);
            auto* mh = new TH1D(("mc_" + histTag + "_" + var).c_str(), "",
                                hc.nbins, hc.xlow, hc.xhigh);
            dh->SetDirectory(nullptr);
            mh->SetDirectory(nullptr);
            dh->SetMarkerStyle(20);
            dh->SetMarkerSize(0.8);
            mh->SetMarkerStyle(21);
            mh->SetMarkerSize(0.8);
            hists.data[var] = dh;
            hists.mc[var] = mh;
        }
    }

    auto fillTree = [&](TTree* tree, bool isData) {
        if (!tree) return;

        BranchBinder b;
        b.bind(tree, ch);
        const Long64_t n = tree->GetEntries();

        for (Long64_t i = 0; i < n; ++i) {
            tree->GetEntry(i);

            if (b.has_runnum && is_excluded_run(b.runnum)) continue;
            if (!(b.has_detector1 && b.has_detector2)) continue;

            bool validTopology = false;
            const Topology topo = topologyFromDetectors(b.detector1, b.detector2, validTopology);
            if (!validTopology) continue;

            const auto cfgIt = globalConfigs.find(topo);
            if (cfgIt == globalConfigs.end()) continue;
            if (!passesGlobalForStage(b, period_label, cfgIt->second)) continue;

            const auto vals = b.valuesMap(ch);
            const auto cumulativeIt = cumulative.find(topo);
            if (cumulativeIt == cumulative.end()) continue;

            const auto& cuts = isData ? cumulativeIt->second.data : cumulativeIt->second.mc;
            if (!passes3SigmaCuts(cuts, vals)) continue;

            FilledHists& target = out.at(topo);
            auto& histMap = isData ? target.data : target.mc;
            for (const auto& kv : cfg) {
                const auto valueIt = vals.find(kv.first);
                if (valueIt != vals.end()) {
                    histMap.at(kv.first)->Fill(valueIt->second);
                }
            }
        }
    };

    fillTree(dataTree, true);
    fillTree(mcTree, false);

    return out;
}

// -------------------- plotting (unchanged aesthetics) --------------------

static std::mutex g_plot_mutex;

static void saveStagePlots(const FilledHists& H, const HistList& cfg, Channel ch,
                           const std::string& prettyPeriod, Topology topo,
                           const std::string& outPlotDir, const std::string& suffix)
{
    std::lock_guard<std::mutex> lock(g_plot_mutex); // serialize ROOT drawing

    int n = static_cast<int>(cfg.size());
    int cols = 4, rows = (n + cols - 1) / cols;
    TCanvas c("c", "", 2400, rows * 640);
    c.Divide(cols, rows, 0.002, 0.002);

    auto isGaussianVar = [&](const std::string& v)->bool {
        return (v == "theta_gamma_gamma" || v == "theta_pi0_pi0" || v == "pTmiss");
    };

    {
        std::string title = channelPretty(ch) + "  |  " + prettyPeriod + "  |  " + topoToKey(topo) + "  |  " + suffix;
        TPaveText* pt = new TPaveText(0.10, 0.955, 0.90, 0.995, "NDC");
        pt->SetFillColor(0);
        pt->SetFillStyle(0);
        pt->SetBorderSize(0);
        pt->SetTextAlign(22);
        pt->SetTextFont(42);
        pt->SetTextSize(0.025);
        pt->AddText(title.c_str());
        pt->Draw();
    }

    int pad = 1;
    for (const auto& kv : cfg) {
        const std::string& var = kv.first;
        TH1D* dh = H.data.at(var);
        TH1D* mh = H.mc.at(var);

        c.cd(pad++);
        gPad->SetLeftMargin(0.21);
        gPad->SetBottomMargin(0.13);
        gPad->SetRightMargin(0.06);
        gPad->SetTopMargin(0.08);
        gPad->SetTickx(1); gPad->SetTicky(1);

        if (dh) { dh->SetLineColor(kBlue); dh->SetMarkerColor(kBlue); dh->SetLineWidth(2); dh->SetStats(0); }
        if (mh) { mh->SetLineColor(kRed);  mh->SetMarkerColor(kRed);  mh->SetLineWidth(2); mh->SetStats(0); }

        normalizeHist(dh); normalizeHist(mh);
        if (dh) {
            dh->GetYaxis()->SetTitle("Normalized counts");
            dh->GetXaxis()->SetTitle(formatLabelName(var, ch).c_str());
            dh->GetXaxis()->SetTitleOffset(1.10);
            dh->GetYaxis()->SetTitleOffset(2.20);
        }
        if (mh && !dh) {
            mh->GetYaxis()->SetTitle("Normalized counts");
            mh->GetXaxis()->SetTitle(formatLabelName(var, ch).c_str());
            mh->GetXaxis()->SetTitleOffset(1.10);
            mh->GetYaxis()->SetTitleOffset(2.20);
        }

        double maxv = 0.0;
        if (dh) maxv = std::max(maxv, dh->GetMaximum());
        if (mh) maxv = std::max(maxv, mh->GetMaximum());
        if (dh) dh->SetMaximum(maxv * 1.25);
        if (!dh && mh) mh->SetMaximum(maxv * 1.25);

        if (dh) dh->Draw("E1");
        if (mh) mh->Draw(dh ? "E1 SAME" : "E1");

        bool isG = (var == "theta_gamma_gamma" || var == "theta_pi0_pi0" || var == "pTmiss");
        double mu_d = 0.0, sg_d = 0.0, mu_m = 0.0, sg_m = 0.0;

        if (isG) {
            if (dh && dh->GetEntries() > 0) {
                auto ms = fitGaussianLeftSide(dh);
                mu_d = ms.first; sg_d = ms.second;
                TF1 f(("fdata_" + var).c_str(), "gaus(0)", dh->GetXaxis()->GetXmin(), dh->GetXaxis()->GetXmax());
                f.SetParameters(dh->GetMaximum(), mu_d, sg_d);
                f.SetLineColor(kBlue + 1); f.SetLineStyle(2); f.SetLineWidth(2); f.Draw("SAME");
            } else { mu_d = dh ? dh->GetMean() : 0.0; sg_d = dh ? dh->GetStdDev() : 0.0; }
            if (mh && mh->GetEntries() > 0) {
                auto ms = fitGaussianLeftSide(mh);
                mu_m = ms.first; sg_m = ms.second;
                TF1 f(("fmc_" + var).c_str(), "gaus(0)", mh->GetXaxis()->GetXmin(), mh->GetXaxis()->GetXmax());
                f.SetParameters(mh->GetMaximum(), mu_m, sg_m);
                f.SetLineColor(kRed + 1); f.SetLineStyle(2); f.SetLineWidth(2); f.Draw("SAME");
            } else { mu_m = mh ? mh->GetMean() : 0.0; sg_m = mh ? mh->GetStdDev() : 0.0; }
        } else {
            mu_d = dh ? dh->GetMean() : 0.0; sg_d = dh ? dh->GetStdDev() : 0.0;
            mu_m = mh ? mh->GetMean() : 0.0; sg_m = mh ? mh->GetStdDev() : 0.0;
        }

        auto dataLine = TString::Format("Data (#mu=%.3f, #sigma=%.3f)", mu_d, sg_d);
        auto mcLine   = TString::Format("MC (#mu=%.3f, #sigma=%.3f)",   mu_m, sg_m);

        TLegend* leg = new TLegend(0.50, 0.68, 0.94, 0.92);
        leg->SetFillStyle(1001);
        leg->SetFillColor(kWhite);
        leg->SetBorderSize(1);
        leg->SetLineColor(kBlack);
        leg->SetTextFont(42);
        leg->SetTextSize(0.030);
        leg->SetMargin(0.12);
        if (dh) leg->AddEntry(dh, dataLine, "lep");
        if (mh) leg->AddEntry(mh, mcLine,   "lep");
        leg->Draw();
    }

    std::string fname = outPlotDir + "/" + prettyPeriod + "_" + topoToKey(topo) + "_" + suffix + "_comparison.png";
    c.SaveAs(fname.c_str());
}


// -------------------- cut statistics --------------------

static double histQuantile(TH1D* h, double q) {
    if (!h || h->GetEntries() <= 0) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    q = std::max(0.0, std::min(1.0, q));
    double prob[1] = {q};
    double quant[1] = {std::numeric_limits<double>::quiet_NaN()};
    h->GetQuantiles(1, quant, prob);
    return quant[0];
}

static Stats makeSymmetricSigmaStats(TH1D* h, double sigmaMultiplier) {
    Stats s = meanStd(h);
    s.mode = "symmetric_3sigma";
    s.quantile = 0.0;

    if (std::isfinite(s.mean) && std::isfinite(s.std) && s.std > 0.0) {
        s.cut_low = s.mean - sigmaMultiplier * s.std;
        s.cut_high = s.mean + sigmaMultiplier * s.std;
    } else {
        s.cut_low = -1.0e99;
        s.cut_high =  1.0e99;
    }

    return s;
}

static Stats makeUpperQuantileStats(TH1D* h, double q) {
    Stats s = meanStd(h);
    s.mode = "upper_quantile";
    s.quantile = q;
    s.cut_low = 0.0;
    s.cut_high = histQuantile(h, q);

    if (!std::isfinite(s.cut_high)) {
        s.cut_high = 1.0e99;
    }

    return s;
}

// -------------------- cumulative updates --------------------

static void updateCumulativeCuts(const FilledHists& H, const StagePlan& stage, CutDict& cumulative, double upperTailQuantile, double symmetricSigmaMultiplier) {
    for (const std::string& var : stage.vars) {
        TH1D* dh = H.data.count(var) ? H.data.at(var) : nullptr;
        TH1D* mh = H.mc.count(var)   ? H.mc.at(var)   : nullptr;

        if (isUpperTailQuantileVar(var)) {
            cumulative.data[var] = makeUpperQuantileStats(dh, upperTailQuantile);
            cumulative.mc[var]   = makeUpperQuantileStats(mh, upperTailQuantile);
        } else {
            cumulative.data[var] = makeSymmetricSigmaStats(dh, symmetricSigmaMultiplier);
            cumulative.mc[var]   = makeSymmetricSigmaStats(mh, symmetricSigmaMultiplier);
        }
    }
}

// -------------------- JSON writers --------------------

static void writeStatsJson(std::ostream& os, const Stats& s) {
    os << "{\"mean\":" << std::setprecision(8) << s.mean
       << ",\"std\":"  << std::setprecision(8) << s.std
       << ",\"cut_low\":" << std::setprecision(8) << s.cut_low
       << ",\"cut_high\":" << std::setprecision(8) << s.cut_high
       << ",\"mode\":\"" << s.mode << "\""
       << ",\"quantile\":" << std::setprecision(8) << s.quantile
       << "}";
}

static void writeCutDictJson(std::ostream& os, const CutDict& cd) {
    os << "{\"data\":{";
    bool first = true;
    for (const auto& kv : cd.data) {
        if (!first) os << ",";
        first = false;
        os << "\"" << kv.first << "\":";
        writeStatsJson(os, kv.second);
    }
    os << "},\"mc\":{";
    first = true;
    for (const auto& kv : cd.mc) {
        if (!first) os << ",";
        first = false;
        os << "\"" << kv.first << "\":";
        writeStatsJson(os, kv.second);
    }
    os << "}}";
}

static double safeRatio(double num, double den) {
    if (!(std::isfinite(num) && std::isfinite(den)) || den == 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    return num / den;
}

static std::string fmtDiag(double v, int precision = 6) {
    if (!std::isfinite(v)) {
        return "nan";
    }

    std::ostringstream ss;
    ss << std::fixed << std::setprecision(precision) << v;
    return ss.str();
}

static const CutDict* findCutDict(const std::map<std::string, CutDict>& combined,
                                  const std::string& prefix,
                                  const std::string& period,
                                  const std::string& topo) {
    const std::string key = prefix + "_" + period + "_" + topo;
    auto it = combined.find(key);

    if (it == combined.end()) {
        return nullptr;
    }

    return &(it->second);
}

static void writeExclusivityRatioDiagnostics(const std::string& outJsonDir,
                                             const std::map<std::string, CutDict>& combined) {
    const std::string path = outJsonDir + "/exclusivity_mc_shape_double_ratios.csv";
    std::ofstream out(path);

    if (!out) {
        std::cerr << "[exclusivity_cuts] WARNING: cannot write diagnostic CSV: " << path << std::endl;
        return;
    }

    out << "channel,topology,variable,"
        << "fa18_inb_mean,sp18_inb_mean,fa18_out_mean,sp18_out_mean,"
        << "sp18_inb_over_fa18_inb_mean,sp18_out_over_fa18_out_mean,mean_double_ratio,"
        << "fa18_inb_std,sp18_inb_std,fa18_out_std,sp18_out_std,"
        << "sp18_inb_over_fa18_inb_std,sp18_out_over_fa18_out_std,std_double_ratio\n";

    const std::vector<std::string> channels = {"DVCS", "eppi0"};
    const std::vector<std::string> topologies = {"FD_FD", "CD_FD", "CD_FT"};

    for (const std::string& ch : channels) {
        for (const std::string& topo : topologies) {
            const CutDict* fa_in = findCutDict(combined, ch, "Fa18_Inb", topo);
            const CutDict* sp_in = findCutDict(combined, ch, "Sp18_Inb", topo);
            const CutDict* fa_out = findCutDict(combined, ch, "Fa18_Out", topo);
            const CutDict* sp_out = findCutDict(combined, ch, "Sp18_Out", topo);

            if (!(fa_in && sp_in && fa_out && sp_out)) {
                continue;
            }

            std::vector<std::string> vars;

            for (const auto& kv : fa_in->mc) {
                vars.push_back(kv.first);
            }

            std::sort(vars.begin(), vars.end());

            for (const std::string& var : vars) {
                auto getStats = [&](const CutDict* cd)->Stats {
                    auto it = cd->mc.find(var);
                    if (it == cd->mc.end()) {
                        return {std::numeric_limits<double>::quiet_NaN(),
                                std::numeric_limits<double>::quiet_NaN()};
                    }
                    return it->second;
                };

                const Stats fi = getStats(fa_in);
                const Stats si = getStats(sp_in);
                const Stats fo = getStats(fa_out);
                const Stats so = getStats(sp_out);

                const double rin_mean = safeRatio(si.mean, fi.mean);
                const double rout_mean = safeRatio(so.mean, fo.mean);
                const double d_mean = safeRatio(rin_mean, rout_mean);

                const double rin_std = safeRatio(si.std, fi.std);
                const double rout_std = safeRatio(so.std, fo.std);
                const double d_std = safeRatio(rin_std, rout_std);

                out << ch << "," << topo << "," << var << ","
                    << fmtDiag(fi.mean) << "," << fmtDiag(si.mean) << ","
                    << fmtDiag(fo.mean) << "," << fmtDiag(so.mean) << ","
                    << fmtDiag(rin_mean) << "," << fmtDiag(rout_mean) << ","
                    << fmtDiag(d_mean) << ","
                    << fmtDiag(fi.std) << "," << fmtDiag(si.std) << ","
                    << fmtDiag(fo.std) << "," << fmtDiag(so.std) << ","
                    << fmtDiag(rin_std) << "," << fmtDiag(rout_std) << ","
                    << fmtDiag(d_std) << "\n";

                if (ch == "DVCS") {
                    std::cout << "[exclusivity_cuts][MC-SHAPE-RATIO] channel=" << ch
                              << " topo=" << topo
                              << " var=" << var
                              << " mean_Sp18Inb/Fa18Inb=" << fmtDiag(rin_mean)
                              << " mean_Sp18Out/Fa18Out=" << fmtDiag(rout_mean)
                              << " mean_double=" << fmtDiag(d_mean)
                              << " std_Sp18Inb/Fa18Inb=" << fmtDiag(rin_std)
                              << " std_Sp18Out/Fa18Out=" << fmtDiag(rout_std)
                              << " std_double=" << fmtDiag(d_std)
                              << std::endl;
                }
            }
        }
    }

    std::cout << "[exclusivity_cuts] Wrote MC shape-ratio diagnostic CSV: " << path << std::endl;
}

static void writeCombinedJson(const std::string& outJsonDir,
                              const std::map<std::string, CutDict>& combined) {
    std::string path = outJsonDir + "/combined_cuts.json";
    std::ofstream ofs(path);
    if (!ofs) {
        std::cerr << "[Error] Cannot open combined JSON: " << path << std::endl;
        return;
    }

    ofs << "{";
    bool firstKey = true;
    for (const auto& kv : combined) {
        if (!firstKey) ofs << ",";
        firstKey = false;
        ofs << "\n  \"" << kv.first << "\": ";
        writeCutDictJson(ofs, kv.second);
    }
    ofs << "\n}\n";
    std::cout << "[Wrote combined JSON] " << path << std::endl;
}

// -------------------- per-period worker --------------------

struct PeriodWork {
    std::string label;  // e.g., "fa18_inb"
    TTree* dvcs_data = nullptr;
    TTree* dvcs_mc   = nullptr;
    TTree* eppi0_data = nullptr;
    TTree* eppi0_mc   = nullptr;
};

static void processOneChannelAllTopologies(
    const std::string& prettyPeriod, Channel ch,
    const std::string& period_label,
    TTree* dataTree, TTree* mcTree,
    const std::string& outPlotDir,
    std::map<Topology, CutDict>& outCuts,
    double upperTailQuantile,
    double symmetricSigmaMultiplier,
    bool makeCutExtractionComparisonPlots)
{
    static const Topology topologies[] = {
        Topology::FD_FD, Topology::CD_FD, Topology::CD_FT
    };

    const auto stages = buildStages(ch);
    const auto cfg = getHistConfigs(ch);

    std::map<Topology, CutDict> cumulative;
    for (Topology topo : topologies) {
        cumulative.emplace(topo, CutDict{});
    }

    const int numStages = static_cast<int>(stages.size()) + 1;
    for (int s = 0; s < numStages; ++s) {
        auto allHists = fillStageHistsAllTopologies(
            dataTree, mcTree, ch, period_label, cumulative, cfg, s);

        for (Topology topo : topologies) {
            FilledHists& H = allHists.at(topo);
            if (makeCutExtractionComparisonPlots) {
                saveStagePlots(H, cfg, ch, prettyPeriod, topo,
                               outPlotDir, "cut_" + std::to_string(s));
            }

            if (s < static_cast<int>(stages.size())) {
                updateCumulativeCuts(H, stages[s], cumulative.at(topo), upperTailQuantile, symmetricSigmaMultiplier);
            }

            for (auto& kv : H.data) delete kv.second;
            for (auto& kv : H.mc) delete kv.second;
        }
    }

    outCuts = std::move(cumulative);
}

static void processPeriod(
    const PeriodWork& W,
    const std::string& outJsonDir,
    const std::string& outPlotDir,
    std::map<std::string, CutDict>& combined_out,
    std::mutex& combined_mutex,
    double upperTailQuantile,
    double symmetricSigmaMultiplier,
    bool makeCutExtractionComparisonPlots)
{
    TH1::AddDirectory(kFALSE);

    if (W.dvcs_data && W.dvcs_mc) {
        const std::string pretty = periodCode(Channel::DVCS, W.label);
        std::map<Topology, CutDict> cutsByTopology;
        processOneChannelAllTopologies(pretty, Channel::DVCS, W.label,
                                       W.dvcs_data, W.dvcs_mc, outPlotDir,
                                       cutsByTopology, upperTailQuantile, symmetricSigmaMultiplier,
                                       makeCutExtractionComparisonPlots);

        {
            std::lock_guard<std::mutex> lock(combined_mutex);
            for (auto& kv : cutsByTopology) {
                combined_out[pretty + "_" + topoToKey(kv.first)] = std::move(kv.second);
            }
        }
        std::cout << "[Done] Exclusivity cuts for " << pretty << std::endl;
    }

    if (W.eppi0_data && W.eppi0_mc) {
        const std::string pretty = periodCode(Channel::EPPI0, W.label);
        std::map<Topology, CutDict> cutsByTopology;
        processOneChannelAllTopologies(pretty, Channel::EPPI0, W.label,
                                       W.eppi0_data, W.eppi0_mc, outPlotDir,
                                       cutsByTopology, upperTailQuantile, symmetricSigmaMultiplier,
                                       makeCutExtractionComparisonPlots);

        {
            std::lock_guard<std::mutex> lock(combined_mutex);
            for (auto& kv : cutsByTopology) {
                combined_out[pretty + "_" + topoToKey(kv.first)] = std::move(kv.second);
            }
        }
        std::cout << "[Done] Exclusivity cuts for " << pretty << std::endl;
    }
}


// -------------------- optional deep diagnostics --------------------

static std::string periodTitleFromCode(const std::string& code) {
    if (code.find("Fa18_Inb") != std::string::npos) return "Fa18 Inb";
    if (code.find("Sp18_Inb") != std::string::npos) return "Sp18 Inb";
    if (code.find("Fa18_Out") != std::string::npos) return "Fa18 Out";
    if (code.find("Sp18_Out") != std::string::npos) return "Sp18 Out";
    if (code.find("Sp19_Inb") != std::string::npos) return "Sp19 Inb";
    return code;
}

static int periodColor(const std::string& code) {
    if (code.find("Fa18_Inb") != std::string::npos) return kBlue + 1;
    if (code.find("Sp18_Inb") != std::string::npos) return kRed + 1;
    if (code.find("Fa18_Out") != std::string::npos) return kGreen + 2;
    if (code.find("Sp18_Out") != std::string::npos) return kMagenta + 1;
    if (code.find("Sp19_Inb") != std::string::npos) return kOrange + 7;
    return kBlack;
}

static int periodMarker(const std::string& code) {
    if (code.find("Fa18_Inb") != std::string::npos) return 20;
    if (code.find("Sp18_Inb") != std::string::npos) return 21;
    if (code.find("Fa18_Out") != std::string::npos) return 22;
    if (code.find("Sp18_Out") != std::string::npos) return 23;
    if (code.find("Sp19_Inb") != std::string::npos) return 33;
    return 20;
}

static GlobalCutConfig globalConfigForTopo(Topology topo, bool bypassGlobalPTmiss) {
    GlobalCutConfig gcfg = default_global_cuts();
    gcfg.enable_topology_filter = true;

    if (topo == Topology::FD_FD) { gcfg.required_detector1 = 1; gcfg.required_detector2 = 1; }
    if (topo == Topology::CD_FD) { gcfg.required_detector1 = 2; gcfg.required_detector2 = 1; }
    if (topo == Topology::CD_FT) { gcfg.required_detector1 = 2; gcfg.required_detector2 = 0; }

    if (bypassGlobalPTmiss) {
        gcfg.enable_pTmiss_cut = false;
        gcfg.pTmiss_max = 1.0e9;
    }

    return gcfg;
}

static bool passesGlobalWithConfig(const BranchBinder& b,
                                   const std::string& period_label,
                                   const GlobalCutConfig& gcfg) {
    if (!(b.has_t1 && b.has_open_angle_ep2)) return false;
    if (gcfg.enable_pTmiss_cut && !b.has_pTmiss) return false;
    if (!(b.has_detector1 && b.has_detector2)) return false;

    if (global_cuts_require_sector_phi(gcfg)) {
        if (!(b.has_e_phi && b.has_p1_phi && b.has_p2_phi)) {
            throw std::runtime_error("[exclusivity_cuts] FATAL: sector selection requires e_phi, p1_phi, p2_phi.");
        }
    }

    if (global_cuts_require_auxiliary_kinematics(gcfg)) {
        if (!(b.has_e_theta && b.has_e_phi &&
              b.has_p1_theta && b.has_p1_phi &&
              b.has_p2_p && b.has_p2_theta && b.has_p2_phi)) {
            throw std::runtime_error("[exclusivity_cuts] FATAL: auxiliary fiducial cuts require e_theta, e_phi, p1_theta, p1_phi, p2_p, p2_theta, p2_phi.");
        }

        return passes_global_cuts(b.t1, b.open_angle_ep2, b.pTmiss,
                                  b.detector1, b.detector2,
                                  period_label,
                                  b.e_p, b.e_theta, b.e_phi,
                                  b.p1_theta, b.p1_phi,
                                  b.p2_p, b.p2_theta, b.p2_phi,
                                  gcfg);
    }

    if (gcfg.enable_dvcsgen_ycol_cut) {
        if (!(b.has_e_p && b.has_e_theta && b.has_e_phi &&
              b.has_p2_p && b.has_p2_theta && b.has_p2_phi)) {
            throw std::runtime_error("[exclusivity_cuts] FATAL: dvcsgen ycol cut requires e_p, e_theta, e_phi, p2_p, p2_theta, p2_phi.");
        }

        if (global_cuts_require_sector_phi(gcfg)) {
            return passes_global_cuts(b.t1, b.open_angle_ep2, b.pTmiss,
                                      b.detector1, b.detector2,
                                      period_label,
                                      b.e_p, b.e_theta, b.e_phi, b.p1_theta, b.p1_phi,
                                      b.p2_p, b.p2_theta, b.p2_phi,
                                      gcfg);
        }

        return passes_global_cuts(b.t1, b.open_angle_ep2, b.pTmiss,
                                  b.detector1, b.detector2,
                                  period_label,
                                  b.e_p, b.e_theta, b.e_phi,
                                  b.p2_p, b.p2_theta, b.p2_phi,
                                  gcfg);
    }

    if (global_cuts_require_sector_phi(gcfg)) {
        return passes_global_cuts(b.t1, b.open_angle_ep2, b.pTmiss,
                                  b.detector1, b.detector2,
                                  b.e_phi, b.p1_phi, b.p2_phi,
                                  gcfg);
    }

    return passes_global_cuts(b.t1, b.open_angle_ep2, b.pTmiss,
                              b.detector1, b.detector2,
                              gcfg);
}

static bool passesCutsExceptVariable(const std::map<std::string, Stats>& cuts,
                                     const std::map<std::string, double>& values,
                                     const std::string& exceptVar) {
    for (const auto& kv : cuts) {
        if (kv.first == exceptVar) continue;

        auto it = values.find(kv.first);
        if (it == values.end()) continue;

        const Stats& s = kv.second;
        if (!(std::isfinite(s.mean) && std::isfinite(s.std) && s.std > 0.0)) continue;

        if (!withinCutWindow(it->second, s)) return false;
    }

    return true;
}

static std::vector<std::string> defaultDiagnosticVariables(Channel ch,
                                                           const ExclusivityDiagnosticConfig& diagCfg) {
    if (!diagCfg.variables.empty()) {
        std::vector<std::string> vars;
        const HistList cfg = getHistConfigs(ch);
        std::set<std::string> allowed;
        for (const auto& kv : cfg) allowed.insert(kv.first);

        for (const std::string& v : diagCfg.variables) {
            if (allowed.count(v)) vars.push_back(v);
        }
        return vars;
    }

    if (ch == Channel::DVCS) {
        return {"pTmiss", "theta_gamma_gamma", "Emiss2", "Mx2_1"};
    }

    return {"pTmiss", "theta_pi0_pi0", "Emiss2", "Mx2_1"};
}

static const HistCfg* findHistCfg(const HistList& cfg, const std::string& var) {
    for (const auto& kv : cfg) {
        if (kv.first == var) return &(kv.second);
    }
    return nullptr;
}

struct DiagnosticHistResult {
    std::unique_ptr<TH1D> hist;
    double selected = 0.0;
    double passSymmetric = 0.0;
    double passOneSidedUpper = 0.0;
    double mean = std::numeric_limits<double>::quiet_NaN();
    double std = std::numeric_limits<double>::quiet_NaN();
};

static DiagnosticHistResult fillMcDiagnosticHist(TTree* tree,
                                                 Channel ch,
                                                 Topology topo,
                                                 const std::string& period_label,
                                                 const std::string& period_code,
                                                 const std::string& var,
                                                 const HistCfg& hc,
                                                 const CutDict& finalCuts,
                                                 bool bypassGlobalPTmiss) {
    DiagnosticHistResult R;
    const std::string hname = "diag_mc_" + channelToStr(ch) + "_" + period_code + "_" + topoToKey(topo) + "_" + var +
                              (bypassGlobalPTmiss ? "_bypassGlobalPT" : "_afterGlobal");
    R.hist.reset(new TH1D(hname.c_str(), "", hc.nbins, hc.xlow, hc.xhigh));
    R.hist->SetDirectory(nullptr);
    R.hist->Sumw2();

    if (!tree) return R;

    GlobalCutConfig gcfg = globalConfigForTopo(topo, bypassGlobalPTmiss);

    BranchBinder b;
    b.bind(tree, ch);

    Long64_t n = tree->GetEntries();
    for (Long64_t i = 0; i < n; ++i) {
        tree->GetEntry(i);

        if (b.has_runnum && is_excluded_run(b.runnum)) continue;
        if (!(b.has_detector1 && b.has_detector2)) continue;
        if (!passesTopology(b.detector1, b.detector2, topo)) continue;
        if (!passesGlobalWithConfig(b, period_label, gcfg)) continue;

        auto vals = b.valuesMap(ch);
        auto it = vals.find(var);
        if (it == vals.end()) continue;
        if (!std::isfinite(it->second)) continue;

        if (!passesCutsExceptVariable(finalCuts.mc, vals, var)) continue;

        R.hist->Fill(it->second);
        R.selected += 1.0;

        auto cs = finalCuts.mc.find(var);
        if (cs != finalCuts.mc.end()) {
            if (withinCutWindow(it->second, cs->second)) R.passSymmetric += 1.0;
            if (it->second <= cs->second.cut_high) R.passOneSidedUpper += 1.0;
        }
    }

    if (R.hist && R.hist->GetEntries() > 0) {
        R.mean = R.hist->GetMean();
        R.std = R.hist->GetStdDev();
    }

    return R;
}

static void drawVerticalLine(double x, double yMax, int color, int style, int width = 2) {
    TLine* line = new TLine(x, 0.0, x, yMax);
    line->SetLineColor(color);
    line->SetLineStyle(style);
    line->SetLineWidth(width);
    line->Draw("SAME");
}

static void runExclusivityOverlayDiagnosticsForChannel(
    Channel ch,
    const std::vector<PeriodWork>& work,
    const std::string& outJsonDir,
    const std::string& outPlotDir,
    const std::map<std::string, CutDict>& combined,
    const ExclusivityDiagnosticConfig& diagCfg)
{
    if (!diagCfg.make_mc_period_overlay_plots && !diagCfg.write_mc_period_overlay_csv) return;

    const std::string diagPlotDir = outPlotDir + "/diagnostics";
    gSystem->mkdir(diagPlotDir.c_str(), true);

    const std::string csvPath = outJsonDir + "/exclusivity_mc_overlay_diagnostic_counts.csv";
    static std::mutex csv_mutex;

    std::ofstream csv;
    if (diagCfg.write_mc_period_overlay_csv) {
        std::lock_guard<std::mutex> lock(csv_mutex);
        const bool exists = !std::ifstream(csvPath).fail();
        csv.open(csvPath, std::ios::app);
        if (!exists) {
            csv << "channel,topology,variable,mode,period,selected,pass_nominal,pass_upper_bound,"
                << "frac_nominal,frac_upper_bound,mean,std,cut_mean,cut_std,cut_low,cut_high,cut_mode,quantile,global_pTmiss_enabled,global_pTmiss_max\n";
        }
    }

    const HistList cfg = getHistConfigs(ch);
    const std::vector<std::string> vars = defaultDiagnosticVariables(ch, diagCfg);
    const std::vector<Topology> topologies = {Topology::FD_FD, Topology::CD_FD, Topology::CD_FT};

    for (Topology topo : topologies) {
        for (const std::string& var : vars) {
            const HistCfg* hc = findHistCfg(cfg, var);
            if (!hc) continue;

            std::vector<bool> bypassModes = {false};
            if (var == "pTmiss" && diagCfg.make_pTmiss_before_global_pTmiss_plots) {
                bypassModes.push_back(true);
            }

            for (bool bypassGlobalPTmiss : bypassModes) {
                const std::string mode = bypassGlobalPTmiss ? "before_global_pTmiss" : "after_global_cuts";

                std::vector<std::string> codes;
                std::vector<DiagnosticHistResult> results;
                std::vector<Stats> cutStats;

                for (const auto& W : work) {
                    TTree* tree = (ch == Channel::DVCS) ? W.dvcs_mc : W.eppi0_mc;
                    if (!tree) continue;

                    const std::string code = periodCode(ch, W.label);
                    const std::string key = code + "_" + topoToKey(topo);
                    auto itCut = combined.find(key);
                    if (itCut == combined.end()) continue;

                    auto itVarCut = itCut->second.mc.find(var);
                    if (itVarCut == itCut->second.mc.end()) continue;

                    DiagnosticHistResult R = fillMcDiagnosticHist(tree, ch, topo, W.label, code, var, *hc, itCut->second, bypassGlobalPTmiss);

                    if (diagCfg.write_mc_period_overlay_csv && csv) {
                        const double fracSym = safeRatio(R.passSymmetric, R.selected);
                        const double fracOne = safeRatio(R.passOneSidedUpper, R.selected);
                        const double lo = itVarCut->second.cut_low;
                        const double hi = itVarCut->second.cut_high;
                        csv << channelPretty(ch) << "," << topoToKey(topo) << "," << var << "," << mode << "," << periodTitleFromCode(code) << ","
                            << fmtDiag(R.selected, 0) << "," << fmtDiag(R.passSymmetric, 0) << "," << fmtDiag(R.passOneSidedUpper, 0) << ","
                            << fmtDiag(fracSym) << "," << fmtDiag(fracOne) << ","
                            << fmtDiag(R.mean) << "," << fmtDiag(R.std) << ","
                            << fmtDiag(itVarCut->second.mean) << "," << fmtDiag(itVarCut->second.std) << ","
                            << fmtDiag(lo) << "," << fmtDiag(hi) << ","
                            << itVarCut->second.mode << "," << fmtDiag(itVarCut->second.quantile) << ","
                            << (default_global_cuts().enable_pTmiss_cut ? "true" : "false") << ","
                            << fmtDiag(default_global_cuts().pTmiss_max) << "\n";
                    }

                    codes.push_back(code);
                    cutStats.push_back(itVarCut->second);
                    results.push_back(std::move(R));
                }

                if (!diagCfg.make_mc_period_overlay_plots || results.empty()) continue;

                std::lock_guard<std::mutex> lock(g_plot_mutex);

                TCanvas c("c_diag", "", 1400, 950);
                c.SetLeftMargin(0.12);
                c.SetBottomMargin(0.12);
                c.SetRightMargin(0.04);
                c.SetTopMargin(0.08);
                c.SetTickx(1);
                c.SetTicky(1);

                double maxv = 0.0;
                for (auto& R : results) {
                    if (!R.hist) continue;
                    normalizeHist(R.hist.get());
                    maxv = std::max(maxv, R.hist->GetMaximum());
                }
                if (maxv <= 0.0) maxv = 1.0;

                bool drawn = false;
                TLegend leg(0.55, 0.60, 0.94, 0.90);
                leg.SetFillColor(kWhite);
                leg.SetBorderSize(1);
                leg.SetTextFont(42);
                leg.SetTextSize(0.030);

                for (size_t i = 0; i < results.size(); ++i) {
                    TH1D* h = results[i].hist.get();
                    if (!h) continue;

                    const int col = periodColor(codes[i]);
                    h->SetLineColor(col);
                    h->SetMarkerColor(col);
                    h->SetMarkerStyle(periodMarker(codes[i]));
                    h->SetMarkerSize(0.7);
                    h->SetLineWidth(2);
                    h->SetStats(0);
                    h->SetMaximum(maxv * 1.45);
                    h->GetYaxis()->SetTitle("Normalized counts");
                    h->GetXaxis()->SetTitle(formatLabelName(var, ch).c_str());
                    h->GetXaxis()->SetTitleOffset(1.05);
                    h->GetYaxis()->SetTitleOffset(1.25);

                    const std::string label = periodTitleFromCode(codes[i]) +
                        TString::Format(" (#mu=%.3g, #sigma=%.3g)", results[i].mean, results[i].std).Data();

                    h->Draw(drawn ? "E1 SAME" : "E1");
                    drawn = true;
                    leg.AddEntry(h, label.c_str(), "lep");
                }

                const double yMax = maxv * 1.35;
                for (size_t i = 0; i < cutStats.size(); ++i) {
                    const int col = periodColor(codes[i]);
                    const double lo = cutStats[i].cut_low;
                    const double hi = cutStats[i].cut_high;

                    if (cutStats[i].mode == "symmetric_3sigma" && diagCfg.draw_symmetric_3sigma_windows) {
                        drawVerticalLine(lo, yMax, col, 3, 1);
                    }
                    if (diagCfg.draw_symmetric_3sigma_windows || diagCfg.draw_one_sided_upper_windows) {
                        drawVerticalLine(hi, yMax, col, 2, 2);
                    }
                }

                if (var == "pTmiss" && diagCfg.draw_global_pTmiss_cut && default_global_cuts().enable_pTmiss_cut) {
                    drawVerticalLine(default_global_cuts().pTmiss_max, yMax, kBlack, 1, 3);
                    leg.AddEntry((TObject*)nullptr, TString::Format("global p_{T}^{miss} < %.3g", default_global_cuts().pTmiss_max), "");
                }

                leg.AddEntry((TObject*)nullptr, "dashed upper lines: active upper bound", "");
                if (diagCfg.draw_symmetric_3sigma_windows) {
                    leg.AddEntry((TObject*)nullptr, "dotted lower lines: symmetric-cut lower bound", "");
                }
                leg.Draw();

                TPaveText title(0.10, 0.925, 0.90, 0.985, "NDC");
                title.SetFillColor(0);
                title.SetFillStyle(0);
                title.SetBorderSize(0);
                title.SetTextAlign(22);
                title.SetTextFont(42);
                title.SetTextSize(0.030);
                const std::string titleText = channelPretty(ch) + " MC diagnostic | " + topoToKey(topo) + " | " + var + " | " + mode;
                title.AddText(titleText.c_str());
                title.Draw();

                const std::string fname = diagPlotDir + "/mc_overlay_" + channelToStr(ch) + "_" + topoToKey(topo) + "_" + var + "_" + mode + ".png";
                c.SaveAs(fname.c_str());

                std::cout << "[exclusivity_cuts][DIAG-PLOT] " << fname << std::endl;
            }
        }
    }

    if (csv) {
        std::cout << "[exclusivity_cuts] Wrote MC overlay diagnostic CSV: " << csvPath << std::endl;
    }
}

static void runOptionalExclusivityDiagnostics(
    const std::vector<PeriodWork>& work,
    const std::string& outJsonDir,
    const std::string& outPlotDir,
    const std::map<std::string, CutDict>& combined,
    const ExclusivityDiagnosticConfig& diagCfg)
{
    if (!diagCfg.enable) return;

    std::cout << "[exclusivity_cuts] Optional deep diagnostics enabled." << std::endl;

    if (diagCfg.include_dvcs) {
        runExclusivityOverlayDiagnosticsForChannel(Channel::DVCS, work, outJsonDir, outPlotDir, combined, diagCfg);
    }

    if (diagCfg.include_eppi0) {
        runExclusivityOverlayDiagnosticsForChannel(Channel::EPPI0, work, outJsonDir, outPlotDir, combined, diagCfg);
    }
}

// -------------------- dispatcher --------------------

void runAllExclusivityCuts(
    const std::map<std::string, TTree*>& dvcsDataTrees,
    const std::map<std::string, TTree*>& dvcsRecMcTrees,
    const std::map<std::string, TTree*>& eppi0DataTrees,
    const std::map<std::string, TTree*>& eppi0RecMcTrees,
    const std::string& outJsonDir,
    const std::string& outPlotDir,
    int maxThreads,
    const ExclusivityDiagnosticConfig& diagCfg)
{
    ROOT::EnableThreadSafety();
    TH1::AddDirectory(kFALSE);
    gStyle->SetOptStat(0);

    if (!outJsonDir.empty()) gSystem->mkdir(outJsonDir.c_str(), true);
    if (!outPlotDir.empty()) gSystem->mkdir(outPlotDir.c_str(), true);

    // Build workers per period
    std::vector<PeriodWork> work;
    work.reserve(CANONICAL_PERIODS().size());

    auto getOrNull = [](const auto& m, const std::string& k)->TTree* {
        auto it = m.find(k);
        return (it != m.end() ? it->second : nullptr);
    };

    auto eppi0CleanBaseFromDvcsBase = [](const std::string& dvcs_base)->std::string {
        // load_trees.cpp now uses explicit EPPI0_* keys for the eppi0 channel:
        //   DVCS_Fa18_inb  -> EPPI0_Fa18_inb_data / EPPI0_Fa18_inb_rec
        // Keep the old DVCS_*_eppi0 lookup below as a compatibility fallback.
        const std::string prefix = "DVCS_";
        if (dvcs_base.rfind(prefix, 0) == 0) {
            return std::string("EPPI0_") + dvcs_base.substr(prefix.size());
        }
        return std::string("EPPI0_") + dvcs_base;
    };

    for (const auto& P : CANONICAL_PERIODS()) {
        const std::string base = P.tree_key;
        const std::string lbl  = P.label;
        const std::string eppi0_base = eppi0CleanBaseFromDvcsBase(base);

        PeriodWork W;
        W.label      = lbl;
        W.dvcs_data  = getOrNull(dvcsDataTrees, base);
        W.dvcs_mc    = getOrNull(dvcsRecMcTrees, base + std::string("_rec"));

        // New clean-tag convention, with old-tag fallbacks for compatibility.
        W.eppi0_data = getOrNull(eppi0DataTrees, eppi0_base + std::string("_data"));
        if (!W.eppi0_data) {
            W.eppi0_data = getOrNull(eppi0DataTrees, base + std::string(SUF_EPPI0));
        }

        W.eppi0_mc = getOrNull(eppi0RecMcTrees, eppi0_base + std::string("_rec"));
        if (!W.eppi0_mc) {
            W.eppi0_mc = getOrNull(eppi0RecMcTrees, base + std::string(SUF_REC_MC));
        }

        if ((W.dvcs_data && W.dvcs_mc) || (W.eppi0_data && W.eppi0_mc)) {
            work.push_back(W);
        }
    }

    if (work.empty()) {
        std::cout << "[Info] No exclusivity jobs found.\n";
        return;
    }

    // Cap to at most 7 threads, and not more than number of periods.
    int nworkers = std::max(1, std::min<int>(maxThreads, 7));
    nworkers = std::min<int>(nworkers, static_cast<int>(work.size()));

    std::map<std::string, CutDict> combined;
    std::mutex combined_mutex;

    std::vector<std::thread> threads;
    threads.reserve(nworkers);

    std::atomic<size_t> idx{0};
    auto worker = [&]() {
        while (true) {
            size_t i = idx.fetch_add(1);
            if (i >= work.size()) break;
            processPeriod(work[i], outJsonDir, outPlotDir, combined, combined_mutex,
                          diagCfg.upper_tail_quantile, diagCfg.symmetric_sigma_multiplier,
                          diagCfg.make_cut_extraction_comparison_plots);
        }
    };

    for (int i = 0; i < nworkers; ++i) threads.emplace_back(worker);
    for (auto& t : threads) t.join();

    writeCombinedJson(outJsonDir, combined);
    writeExclusivityRatioDiagnostics(outJsonDir, combined);
    runOptionalExclusivityDiagnostics(work, outJsonDir, outPlotDir, combined, diagCfg);

    std::cout << "[All done] Exclusivity cuts ran for " << work.size()
              << " period(s) with up to " << nworkers << " worker(s)." << std::endl;
}