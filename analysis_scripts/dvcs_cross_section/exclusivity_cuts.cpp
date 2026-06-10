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
#include <stdexcept>
#include <string>
#include <thread>
#include <chrono>
#include <utility>
#include <vector>

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
    if (ch == Channel::DVCS) stages.push_back(StagePlan{{"pTmiss", "xF", "theta_gamma_gamma"}});
    else                     stages.push_back(StagePlan{{"pTmiss", "xF", "theta_pi0_pi0"}});
    return stages; // final pass is len(stages)+1 in driver
}

// -------------------- topology --------------------

static bool passesTopology(int detector1, int detector2, Topology topo) {
    if (topo == Topology::FD_FD) return (detector1 == 1 && detector2 == 1);
    if (topo == Topology::CD_FD) return (detector1 == 2 && detector2 == 1);
    if (topo == Topology::CD_FT) return (detector1 == 2 && detector2 == 0);
    return false;
}

static bool within3Sigma(double val, const Stats& s) {
    return (val >= s.mean - 3.0 * s.std) && (val <= s.mean + 3.0 * s.std);
}

static bool passes3SigmaCuts(const std::map<std::string, Stats>& cuts,
                             const std::map<std::string, double>& values) {
    for (const auto& kv : cuts) {
        auto it = values.find(kv.first);
        if (it == values.end()) continue;
        if (!within3Sigma(it->second, kv.second)) return false;
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
        bD("p1_phi",  &p1_phi,  has_p1_phi);

        bD("p2_p",     &p2_p,     has_p2_p);
        bD("p2_theta", &p2_theta, has_p2_theta);
        bD("p2_phi",   &p2_phi,   has_p2_phi);
    }

    std::map<std::string, double> valuesMap(Channel ch) const {
        std::map<std::string, double> m;
        if (has_Delta_phi) m["Delta_phi"] = Delta_phi;
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

static FilledHists fillStageHists(
    TTree* dataTree, TTree* mcTree, Topology topo, Channel ch,
    const std::string& period_label,
    const CutDict& cumulative, const HistList& cfg, int stage_index)
{
    FilledHists out;

    const std::string histTag = channelToStr(ch) + "_" + period_label + "_" + topoToKey(topo) + "_stage" + std::to_string(stage_index);

    // Create detached histograms with globally unique names.
    for (const auto& kv : cfg) {
        const std::string& var = kv.first; const HistCfg& hc = kv.second;
        auto* dh = new TH1D(("data_" + histTag + "_" + var).c_str(), "", hc.nbins, hc.xlow, hc.xhigh);
        auto* mh = new TH1D(("mc_"   + histTag + "_" + var).c_str(), "", hc.nbins, hc.xlow, hc.xhigh);
        dh->SetDirectory(nullptr); mh->SetDirectory(nullptr);
        dh->SetMarkerStyle(20); dh->SetMarkerSize(0.8);
        mh->SetMarkerStyle(21); mh->SetMarkerSize(0.8);
        out.data[var] = dh; out.mc[var] = mh;
    }

    // IMPORTANT:
    // This module explicitly loops over Topology and already applies passesTopology().
    // However, default_global_cuts() may have enable_topology_filter=true, and the
    // fail-fast policy requires that we *never* call the legacy overloads in that case.
    //
    // Therefore we:
    //   - make a local mutable cfg copy,
    //   - explicitly set required_detector1/2 to match the current topo,
    //   - and always call the topology-aware passes_global_cuts overloads.
    GlobalCutConfig gcfg = default_global_cuts();
    gcfg.enable_topology_filter = true;

    if (topo == Topology::FD_FD) { gcfg.required_detector1 = 1; gcfg.required_detector2 = 1; }
    if (topo == Topology::CD_FD) { gcfg.required_detector1 = 2; gcfg.required_detector2 = 1; }
    if (topo == Topology::CD_FT) { gcfg.required_detector1 = 2; gcfg.required_detector2 = 0; }

    auto passesGlobal = [&](const BranchBinder& b)->bool {
        if (!(b.has_t1 && b.has_open_angle_ep2 && b.has_pTmiss)) return false;
        if (!(b.has_detector1 && b.has_detector2)) return false;

        if (global_cuts_require_sector_phi(gcfg)) {
            if (!(b.has_e_phi && b.has_p1_phi && b.has_p2_phi)) {
                throw std::runtime_error("[exclusivity_cuts] FATAL: sector selection requires e_phi, p1_phi, p2_phi.");
            }
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
                                          b.e_p, b.e_theta, b.e_phi, b.p1_phi,
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
    };

    // Data loop
    if (dataTree) {
        BranchBinder b; b.bind(dataTree, ch);
        Long64_t n = dataTree->GetEntries();
        for (Long64_t i = 0; i < n; ++i) {
            dataTree->GetEntry(i);

            // Drop excluded runs if runnum is available (GLOBAL blacklist).
            if (b.has_runnum && is_excluded_run(b.runnum)) continue;

            if (!(b.has_detector1 && b.has_detector2)) continue;
            if (!passesTopology(b.detector1, b.detector2, topo)) continue;

            if (!passesGlobal(b)) continue;

            auto vals = b.valuesMap(ch);
            if (!passes3SigmaCuts(cumulative.data, vals)) continue;

            for (const auto& kv : cfg) {
                auto it = vals.find(kv.first);
                if (it != vals.end()) out.data[kv.first]->Fill(it->second);
            }
        }
    }

    // MC loop
    if (mcTree) {
        BranchBinder b; b.bind(mcTree, ch);
        Long64_t n = mcTree->GetEntries();
        for (Long64_t i = 0; i < n; ++i) {
            mcTree->GetEntry(i);

            // Drop excluded runs if runnum is available (usually only for data,
            // but this check is harmless for MC).
            if (b.has_runnum && is_excluded_run(b.runnum)) continue;

            if (!(b.has_detector1 && b.has_detector2)) continue;
            if (!passesTopology(b.detector1, b.detector2, topo)) continue;

            if (!passesGlobal(b)) continue;

            auto vals = b.valuesMap(ch);
            if (!passes3SigmaCuts(cumulative.mc, vals)) continue;

            for (const auto& kv : cfg) {
                auto it = vals.find(kv.first);
                if (it != vals.end()) out.mc[kv.first]->Fill(it->second);
            }
        }
    }

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

// -------------------- cumulative updates --------------------

static void updateCumulativeCuts(const FilledHists& H, const StagePlan& stage, CutDict& cumulative) {
    auto isGaussianVar = [&](const std::string& v)->bool {
        return (v == "theta_gamma_gamma" || v == "theta_pi0_pi0" || v == "pTmiss");
    };
    for (const std::string& var : stage.vars) {
        TH1D* dh = H.data.count(var) ? H.data.at(var) : nullptr;
        TH1D* mh = H.mc.count(var)   ? H.mc.at(var)   : nullptr;
        if (isGaussianVar(var)) {
            auto d = fitGaussianLeftSide(dh);
            auto m = fitGaussianLeftSide(mh);
            cumulative.data[var] = {d.first, std::abs(d.second)};
            cumulative.mc[var]   = {m.first, std::abs(m.second)};
        } else {
            cumulative.data[var] = meanStd(dh);
            cumulative.mc[var]   = meanStd(mh);
        }
    }
}

// -------------------- JSON writers --------------------

static void writeStatsJson(std::ostream& os, const Stats& s) {
    os << "{\"mean\":" << std::setprecision(8) << s.mean
       << ",\"std\":"  << std::setprecision(8) << s.std << "}";
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

static void processOneChannelOneTopology(
    const std::string& prettyPeriod, Channel ch, Topology topo,
    const std::string& period_label,
    TTree* dataTree, TTree* mcTree,
    const std::string& outPlotDir,
    CutDict& outCutsForTopo)
{
    auto stages = buildStages(ch);
    auto cfg    = getHistConfigs(ch);

    CutDict cumulative;
    int numStages = static_cast<int>(stages.size()) + 1;
    for (int s = 0; s < numStages; ++s) {
        auto H = fillStageHists(dataTree, mcTree, topo, ch, period_label, cumulative, cfg, s);
        saveStagePlots(H, cfg, ch, prettyPeriod, topo, outPlotDir, "cut_" + std::to_string(s));
        if (s < static_cast<int>(stages.size())) {
            updateCumulativeCuts(H, stages[s], cumulative);
        }
        for (auto& kv : H.data) delete kv.second;
        for (auto& kv : H.mc)   delete kv.second;
    }

    outCutsForTopo = cumulative;
}

static void processPeriod(
    const PeriodWork& W,
    const std::string& outJsonDir,
    const std::string& outPlotDir,
    std::map<std::string, CutDict>& combined_out,
    std::mutex& combined_mutex)
{
    TH1::AddDirectory(kFALSE);

    if (W.dvcs_data && W.dvcs_mc) {
        std::string pretty = periodCode(Channel::DVCS, W.label);
        for (Topology topo : {Topology::FD_FD, Topology::CD_FD, Topology::CD_FT}) {
            CutDict cutsDVCS;
            processOneChannelOneTopology(pretty, Channel::DVCS, topo, W.label,
                                         W.dvcs_data, W.dvcs_mc, outPlotDir, cutsDVCS);
            std::lock_guard<std::mutex> lock(combined_mutex);
            combined_out[pretty + "_" + topoToKey(topo)] = cutsDVCS;
        }
        std::cout << "[Done] Exclusivity cuts for " << pretty << std::endl;
    }

    if (W.eppi0_data && W.eppi0_mc) {
        std::string pretty = periodCode(Channel::EPPI0, W.label);
        for (Topology topo : {Topology::FD_FD, Topology::CD_FD, Topology::CD_FT}) {
            CutDict cutsPI0;
            processOneChannelOneTopology(pretty, Channel::EPPI0, topo, W.label,
                                         W.eppi0_data, W.eppi0_mc, outPlotDir, cutsPI0);
            std::lock_guard<std::mutex> lock(combined_mutex);
            combined_out[pretty + "_" + topoToKey(topo)] = cutsPI0;
        }
        std::cout << "[Done] Exclusivity cuts for " << pretty << std::endl;
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
    int maxThreads)
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

    // Cap to at most 5 threads, and not more than number of periods.
    int nworkers = std::max(1, std::min<int>(maxThreads, 5));
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
            processPeriod(work[i], outJsonDir, outPlotDir, combined, combined_mutex);
        }
    };

    for (int i = 0; i < nworkers; ++i) threads.emplace_back(worker);
    for (auto& t : threads) t.join();

    writeCombinedJson(outJsonDir, combined);

    std::cout << "[All done] Exclusivity cuts ran for " << work.size()
              << " period(s) with up to " << nworkers << " worker(s)." << std::endl;
}