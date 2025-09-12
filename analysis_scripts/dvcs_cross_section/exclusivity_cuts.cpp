// exclusivity_cuts.cpp
#include "exclusivity_cuts.h"

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

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <thread>
#include <utility>
#include <vector>

// -------------------- helpers: strings and keys --------------------

static std::string channelToStr(Channel ch) { return (ch == Channel::DVCS) ? "dvcs" : "eppi0"; }
static std::string channelPretty(Channel ch) { return ch == Channel::DVCS ? "DVCS" : "eppi0"; }

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
    for (size_t i = 0; i + 1 < nice.size(); ++i) if (nice[i] == '_' && i + 1 < nice.size())
        nice[i + 1] = std::toupper(nice[i + 1]);
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

// -------------------- topology and global kinematic cuts --------------------
//
// Global kinematic cuts (from kin_cuts.py):
// 1) open_angle_ep2 > 5.0 deg
// 2) (-t1) <= 1.0
// 3) pTmiss <= 0.20

static bool passesTopology(int detector1, int detector2, Topology topo) {
    if (topo == Topology::FD_FD) return (detector1 == 1 && detector2 == 1);
    if (topo == Topology::CD_FD) return (detector1 == 2 && detector2 == 1);
    if (topo == Topology::CD_FT) return (detector1 == 2 && detector2 == 0);
    return false;
}

static bool applyKinematicCuts(double t1, double open_angle_ep2, double pTmiss) {
    if (open_angle_ep2 <= 5.0) return false;
    if ((-t1) > 1.0) return false;
    if (pTmiss > 0.20) return false;
    return true;
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

// -------------------- branch binder (exact names) --------------------

struct BranchBinder {
    int detector1 = 0, detector2 = 0; bool has_detector1 = false, has_detector2 = false;

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

    void bind(TTree* t) {
        auto bindI = [&](const char* n, int* a, bool& f){ if (t->GetBranch(n)) { t->SetBranchAddress(n, a); f = true; } };
        auto bindD = [&](const char* n, double* a, bool& f){ if (t->GetBranch(n)) { t->SetBranchAddress(n, a); f = true; } };

        bindI("detector1", &detector1, has_detector1);
        bindI("detector2", &detector2, has_detector2);

        bindD("t1", &t1, has_t1);
        bindD("open_angle_ep2", &open_angle_ep2, has_open_angle_ep2);
        bindD("Emiss2", &Emiss2, has_Emiss2);
        bindD("Mx2", &Mx2, has_Mx2);
        bindD("Mx2_1", &Mx2_1, has_Mx2_1);
        bindD("Mx2_2", &Mx2_2, has_Mx2_2);
        bindD("pTmiss", &pTmiss, has_pTmiss);
        bindD("xF", &xF, has_xF);
        bindD("Delta_phi", &Delta_phi, has_Delta_phi);
        bindD("theta_gamma_gamma", &theta_gamma_gamma, has_theta_gamma_gamma);
        bindD("theta_pi0_pi0", &theta_pi0_pi0, has_theta_pi0_pi0);
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
    double xmin = h->GetXaxis()->GetXmin();
    int peakBin = h->GetMaximumBin();
    double xpeak = h->GetXaxis()->GetBinCenter(peakBin);
    double xr = std::max(xmin + 1e-6, 1.35 * xpeak);

    TF1 f("gausLeft", "gaus(0)", xmin, xr);
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

// -------------------- pretty x-axis labels (ASCII TLatex) --------------------

static std::string formatLabelName(const std::string& var, Channel ch) {
    if (var == "Delta_phi")          return "#Delta#phi (rad)";
    if (var == "theta_gamma_gamma")  return "#theta_{#gamma#gamma} (rad)";
    if (var == "theta_pi0_pi0")      return "#theta_{#pi^{0}#pi^{0}} (rad)";
    if (var == "pTmiss")             return "p_{T}^{miss} (GeV)";
    if (var == "xF")                 return "x_{F}";
    if (var == "Emiss2")             return "E_{miss}^{2} (GeV^{2})";
    if (var == "Mx2")                return "M_{x}^{2} (GeV^{2})";
    if (var == "Mx2_1")              return "M_{x1}^{2} (GeV^{2})";
    if (var == "Mx2_2")              return "M_{x2}^{2} (GeV^{2})";
    return var;
}

// -------------------- stage fill --------------------

struct FilledHists { std::map<std::string, TH1D*> data; std::map<std::string, TH1D*> mc; };

static FilledHists fillStageHists(
    TTree* dataTree, TTree* mcTree, Topology topo, Channel ch,
    const CutDict& cumulative, const HistList& cfg, int stage_index)
{
    FilledHists out;

    // Create hists (detached from any directory to avoid name collisions)
    for (const auto& kv : cfg) {
        const std::string& var = kv.first; const HistCfg& hc = kv.second;
        auto* dh = new TH1D(("data_" + var + "_stage" + std::to_string(stage_index)).c_str(), "", hc.nbins, hc.xlow, hc.xhigh);
        auto* mh = new TH1D(("mc_"   + var + "_stage" + std::to_string(stage_index)).c_str(), "", hc.nbins, hc.xlow, hc.xhigh);
        dh->SetDirectory(nullptr); mh->SetDirectory(nullptr);
        dh->SetMarkerStyle(20); dh->SetMarkerSize(0.8);
        mh->SetMarkerStyle(21); mh->SetMarkerSize(0.8);
        out.data[var] = dh; out.mc[var] = mh;
    }

    // Data loop
    if (dataTree) {
        BranchBinder b; b.bind(dataTree);
        Long64_t n = dataTree->GetEntries();
        for (Long64_t i = 0; i < n; ++i) {
            dataTree->GetEntry(i);
            if (!(b.has_detector1 && b.has_detector2)) continue;
            if (!passesTopology(b.detector1, b.detector2, topo)) continue;
            if (!(b.has_t1 && b.has_open_angle_ep2 && b.has_pTmiss)) continue;
            if (!applyKinematicCuts(b.t1, b.open_angle_ep2, b.pTmiss)) continue;

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
        BranchBinder b; b.bind(mcTree);
        Long64_t n = mcTree->GetEntries();
        for (Long64_t i = 0; i < n; ++i) {
            mcTree->GetEntry(i);
            if (!(b.has_detector1 && b.has_detector2)) continue;
            if (!passesTopology(b.detector1, b.detector2, topo)) continue;
            if (!(b.has_t1 && b.has_open_angle_ep2 && b.has_pTmiss)) continue;
            if (!applyKinematicCuts(b.t1, b.open_angle_ep2, b.pTmiss)) continue;

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

// -------------------- plotting --------------------

static void saveStagePlots(const FilledHists& H, const HistList& cfg, Channel ch,
                           const std::string& prettyPeriod, Topology topo,
                           const std::string& outPlotDir, const std::string& suffix)
{
    // Canvas: grid of pads; add a global title; remove stat box.
    int n = static_cast<int>(cfg.size());
    int cols = 4, rows = (n + cols - 1) / cols;
    TCanvas c("c", "", 2400, rows * 640);
    c.Divide(cols, rows, 0.002, 0.002);

    auto isGaussianVar = [&](const std::string& v)->bool {
        return (v == "theta_gamma_gamma" || v == "theta_pi0_pi0" || v == "pTmiss");
    };

    // Global title (channel, period, topology, stage)
    {
        std::string title = channelPretty(ch) + "  |  " + prettyPeriod + "  |  " + topoToKey(topo) + "  |  " + suffix;
        TPaveText* pt = new TPaveText(0.10, 0.955, 0.90, 0.995, "NDC");
        pt->SetFillColor(0);
        pt->SetFillStyle(0);
        pt->SetBorderSize(0);
        pt->SetTextAlign(22); // center
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
        gPad->SetLeftMargin(0.21);    // more left space to separate Y title
        gPad->SetBottomMargin(0.13);
        gPad->SetRightMargin(0.06);   // a bit of right padding so legend fits
        gPad->SetTopMargin(0.08);
        gPad->SetTickx(1); gPad->SetTicky(1);

        if (dh) { dh->SetLineColor(kBlue); dh->SetMarkerColor(kBlue); dh->SetLineWidth(2); dh->SetStats(0); }
        if (mh) { mh->SetLineColor(kRed);  mh->SetMarkerColor(kRed);  mh->SetLineWidth(2); mh->SetStats(0); }

        // Normalize and set axis titles
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

        // Determine max for nice scaling
        double maxv = 0.0;
        if (dh) maxv = std::max(maxv, dh->GetMaximum());
        if (mh) maxv = std::max(maxv, mh->GetMaximum());
        if (dh) dh->SetMaximum(maxv * 1.25);
        if (!dh && mh) mh->SetMaximum(maxv * 1.25);

        // Draw histograms
        if (dh) dh->Draw("E1");
        if (mh) mh->Draw(dh ? "E1 SAME" : "E1");

        // Compute mu/sigma for legend and optionally draw fits
        double mu_d = 0.0, sg_d = 0.0, mu_m = 0.0, sg_m = 0.0;

        if (isGaussianVar(var)) {
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

        // Legend: boxed, solid white; bottom-left moved left, top-right flush with pad border
        auto dataLine = TString::Format("Data (mu=%.3f, sigma=%.3f)", mu_d, sg_d);
        auto mcLine   = TString::Format("MC (mu=%.3f, sigma=%.3f)",   mu_m, sg_m);

        TLegend* leg = new TLegend(0.50, 0.68, 0.98, 0.88);
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

static void saveFinalCutsJson(const std::string& periodCodeStr, Topology topo,
                              const std::string& outJsonDir, const CutDict& cuts) {
    std::string path = outJsonDir + "/cuts_" + periodCodeStr + "_" + topoToKey(topo) + "_final.json";
    std::ofstream ofs(path);
    if (!ofs) { std::cerr << "[Error] Cannot open JSON: " << path << std::endl; return; }
    writeCutDictJson(ofs, cuts);
    ofs << std::endl;
    std::cout << "[Wrote JSON] " << path << std::endl;
}

static void writeCombinedJson(const std::string& outJsonDir,
                              const std::map<std::string, CutDict>& combined)
{
    std::string path = outJsonDir + "/combined_cuts.json";
    std::ofstream ofs(path);
    if (!ofs) { std::cerr << "[Error] Cannot open combined JSON: " << path << std::endl; return; }

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

// -------------------- per-period driver --------------------

static void runExclusivityCutsSingle(const std::string& runTag, Channel ch,
                                     TTree* dataTree, TTree* mcTree,
                                     const std::string& outJsonDir,
                                     const std::string& outPlotDir,
                                     std::map<std::string, CutDict>& combinedOut)
{
    if (!dataTree || !mcTree) {
        std::cerr << "[Skip] Missing data or MC for " << channelToStr(ch)
                  << " " << runTag << std::endl;
        return;
    }

    std::string pretty = periodCode(ch, runTag);
    auto stages = buildStages(ch);
    auto cfg = getHistConfigs(ch);

    for (Topology topo : {Topology::FD_FD, Topology::CD_FD, Topology::CD_FT}) {
        CutDict cumulative;
        int numStages = static_cast<int>(stages.size()) + 1; // final pass
        for (int s = 0; s < numStages; ++s) {
            auto H = fillStageHists(dataTree, mcTree, topo, ch, cumulative, cfg, s);
            saveStagePlots(H, cfg, ch, pretty, topo, outPlotDir, "cut_" + std::to_string(s));
            if (s < static_cast<int>(stages.size())) updateCumulativeCuts(H, stages[s], cumulative);
            for (auto& kv : H.data) delete kv.second;
            for (auto& kv : H.mc)   delete kv.second;
        }
        saveFinalCutsJson(pretty, topo, outJsonDir, cumulative);
        std::string combinedKey = pretty + "_" + topoToKey(topo);
        combinedOut[combinedKey] = cumulative;
    }

    std::cout << "[Done] Exclusivity cuts for " << pretty << std::endl;
}

// -------------------- dispatcher (single-threaded) --------------------

void runAllExclusivityCuts(
    const std::map<std::string, TTree*>& dvcsDataTrees,
    const std::map<std::string, TTree*>& dvcsRecMcTrees,
    const std::map<std::string, TTree*>& eppi0DataTrees,
    const std::map<std::string, TTree*>& eppi0RecMcTrees,
    const std::string& outJsonDir,
    const std::string& outPlotDir,
    int maxThreads)
{
    TH1::AddDirectory(kFALSE);
    gStyle->SetOptStat(0);

    const std::vector<std::string> runTags = {
        "sp18_inb", "sp18_out", "fa18_inb_supp", "fa18_inb", "fa18_out", "sp19_inb"
    };

    struct Job { std::string runTag; Channel ch; TTree* data; TTree* mc; };
    std::vector<Job> jobs;

    auto getOrNull = [](const auto& m, const std::string& k)->TTree* {
        auto it = m.find(k); return (it != m.end() ? it->second : nullptr);
    };

    // DVCS
    for (const auto& tag : runTags) {
        TTree* data = getOrNull(dvcsDataTrees, tag);
        TTree* mc   = getOrNull(dvcsRecMcTrees, tag + "_rec");
        if (data && mc) jobs.push_back({tag, Channel::DVCS, data, mc});
    }
    // eppi0
    for (const auto& tag : runTags) {
        TTree* data = getOrNull(eppi0DataTrees, tag + "_eppi0");
        TTree* mc   = getOrNull(eppi0RecMcTrees, tag + "_rec_mc");
        if (data && mc) jobs.push_back({tag, Channel::EPPI0, data, mc});
    }

    if (jobs.empty()) { std::cout << "[Info] No exclusivity jobs found.\n"; return; }

    (void)maxThreads; // unused, single-threaded run

    std::map<std::string, CutDict> combined;
    for (const auto& job : jobs) {
        runExclusivityCutsSingle(job.runTag, job.ch, job.data, job.mc, outJsonDir, outPlotDir, combined);
    }
    writeCombinedJson(outJsonDir, combined);

    std::cout << "[All done] Exclusivity cuts ran for " << jobs.size() << " job(s)." << std::endl;
}