// exclusivity_cuts.h
// Multi-stage exclusivity cut extraction with global kinematic cuts.
// - Mirrors your Python flow (stages, topology, data/MC, JSON outputs).
// - Applies global cuts: open_angle_ep2 > 5 deg, (-t) <= 1.0, pTmiss <= 0.20.
// - Writes per-topology JSON summaries and stage plots.
// - Can run up to 3 threads across independent run periods.
//
// You can plug in more detailed pre-cuts later if needed.

#ifndef EXCLUSIVITY_CUTS_H
#define EXCLUSIVITY_CUTS_H

#include <map>
#include <string>
#include <vector>

class TTree;

enum class Channel { DVCS, EPPI0 };
enum class Topology { FD_FD, CD_FD, CD_FT };

struct Stats { double mean = 0.0; double std = 0.0; };

struct CutDict {
    std::map<std::string, Stats> data;
    std::map<std::string, Stats> mc;
};

struct HistCfg { int nbins = 100; double xlow = 0.0; double xhigh = 1.0; };

struct StagePlan { std::vector<std::string> vars; };

// Run multi-stage exclusivity cuts across all periods and both channels.
// Inputs: tree maps produced by loadTrees().
// Outputs: JSON files in outJsonDir, plots in outPlotDir.
// Concurrency: up to maxThreads (capped at 3).
void runAllExclusivityCuts(
    const std::map<std::string, TTree*>& dvcsDataTrees,
    const std::map<std::string, TTree*>& dvcsRecMcTrees,
    const std::map<std::string, TTree*>& eppi0DataTrees,
    const std::map<std::string, TTree*>& eppi0RecMcTrees,
    const std::string& outJsonDir,
    const std::string& outPlotDir,
    int maxThreads = 3);

#endif // EXCLUSIVITY_CUTS_H