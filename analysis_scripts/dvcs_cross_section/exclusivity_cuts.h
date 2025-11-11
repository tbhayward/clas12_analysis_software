// exclusivity_cuts.h
// Multi-stage exclusivity cut extraction with global kinematic cuts and plotting.
// - Universal global cuts (applied to all entries before stage logic):
//     (-t1) < 1.0, open_angle_ep2 > 5.0 deg, pTmiss <= 0.20
// - Stagewise mu/sigma extraction (left-side Gaussian for theta_* and pTmiss)
// - One combined JSON for all periods and topologies: output/jsons/combined_cuts.json
// - Plots are identical to previous version; filenames and styles unchanged.
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

void runAllExclusivityCuts(
    const std::map<std::string, TTree*>& dvcsDataTrees,
    const std::map<std::string, TTree*>& dvcsRecMcTrees,
    const std::map<std::string, TTree*>& eppi0DataTrees,
    const std::map<std::string, TTree*>& eppi0RecMcTrees,
    const std::string& outJsonDir,  // base json dir ("output/jsons")
    const std::string& outPlotDir,
    int maxThreads = 5               // capped internally to at most 5
);

#endif // EXCLUSIVITY_CUTS_H