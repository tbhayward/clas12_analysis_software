// exclusivity_cuts.h
// Multi-stage exclusivity cut extraction with global kinematic cuts.
// - Global cuts: open_angle_ep2 > 5 deg, (-t1) <= 1.0, pTmiss <= 0.20
// - Stagewise mu/sigma extraction (left-side Gaussian for theta_* and pTmiss)
// - JSON outputs per topology, stage comparison plots
// - Threading set to 1 for robust ROOT I/O (we can revisit with RDataFrame)
//
// ASCII-only.

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

void runAllExclusivityCuts(
    const std::map<std::string, TTree*>& dvcsDataTrees,
    const std::map<std::string, TTree*>& dvcsRecMcTrees,
    const std::map<std::string, TTree*>& eppi0DataTrees,
    const std::map<std::string, TTree*>& eppi0RecMcTrees,
    const std::string& outJsonDir,
    const std::string& outPlotDir,
    int maxThreads = 1 // keep 1 for stability; we can revisit
);

#endif // EXCLUSIVITY_CUTS_H