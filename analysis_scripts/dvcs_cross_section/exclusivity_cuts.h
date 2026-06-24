// exclusivity_cuts.h
// Multi-stage exclusivity cut extraction with global kinematic cuts and plotting.
// - Universal global cuts are provided by global_cuts.h.
// - Stagewise mu/sigma extraction.
// - One combined JSON for all periods and topologies: output/jsons/combined_cuts.json.
// - Optional deep diagnostic plots for MC exclusivity variables and cut boundaries.
// - Parallelized by period with a hard cap of 5 threads.

#ifndef EXCLUSIVITY_CUTS_H
#define EXCLUSIVITY_CUTS_H

#include <map>
#include <string>
#include <vector>

class TTree;

enum class Channel { DVCS, EPPI0 };
enum class Topology { FD_FD, CD_FD, CD_FT };

struct Stats { double mean = 0.0; double std = 0.0; };
struct CutDict { std::map<std::string, Stats> data, mc; };
struct HistCfg { int nbins = 100; double xlow = 0.0; double xhigh = 1.0; };
struct StagePlan { std::vector<std::string> vars; };

struct ExclusivityDiagnosticConfig {
    bool enable = false;

    // Recommended default for the present DVCS acceptance investigation.
    bool include_dvcs = true;
    bool include_eppi0 = false;

    // MC-only period overlays are the most direct diagnostic for acceptance.
    bool make_mc_period_overlay_plots = true;
    bool write_mc_period_overlay_csv = true;

    // Draw final period-specific MC cut boundaries on the overlay plots.
    bool draw_symmetric_3sigma_windows = true;
    bool draw_one_sided_upper_windows = true;
    bool draw_global_pTmiss_cut = true;

    // For pTmiss only, also make a version where the early global pTmiss cut is bypassed.
    // This is useful because pTmiss is both a global pre-cut and a later exclusivity cut.
    bool make_pTmiss_before_global_pTmiss_plots = true;

    // If empty, defaults are chosen internally:
    // DVCS: pTmiss, theta_gamma_gamma, Emiss2, Mx2_1
    // eppi0: pTmiss, theta_pi0_pi0, Emiss2, Mx2_1
    std::vector<std::string> variables;
};

void runAllExclusivityCuts(
    const std::map<std::string, TTree*>& dvcsDataTrees,
    const std::map<std::string, TTree*>& dvcsRecMcTrees,
    const std::map<std::string, TTree*>& eppi0DataTrees,
    const std::map<std::string, TTree*>& eppi0RecMcTrees,
    const std::string& outJsonDir,  // base json dir ("output/jsons")
    const std::string& outPlotDir,
    int maxThreads = 5,             // capped internally to at most 5
    const ExclusivityDiagnosticConfig& diagCfg = ExclusivityDiagnosticConfig()
);

#endif // EXCLUSIVITY_CUTS_H
