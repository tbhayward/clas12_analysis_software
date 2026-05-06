#ifndef LOAD_TREES_H
#define LOAD_TREES_H

#include <map>
#include <string>

class TTree;

/**
 * loadTrees
 *
 * Loads all ROOT trees used by the DVCS pass-2 cross-section workflow.
 *
 * Fa18 Inb supplemental data are intentionally excluded from this loader and
 * from the standard analysis workflow.
 *
 * The original maps are preserved:
 *
 *   dataTrees
 *     DVCS data trees.
 *
 *   genMcTrees / recMcTrees
 *     The reference-current non-radiative DVCS MC trees used for the standard
 *     acceptance correction:
 *
 *       Sp18 Inb -> 50 nA
 *       Sp18 Out -> 45 nA
 *       Fa18 Inb -> 50 nA
 *       Fa18 Out -> 50 nA
 *       Sp19 Inb -> 50 nA
 *
 *   eppi0DataTrees / eppi0GenMcTrees / eppi0RecMcTrees
 *     eppi0 data and aaogen MC trees used for the eppi0 normalization study.
 *
 *   eppi0BkgTrees
 *     eppi0 -> DVCS background MC trees used for pi0 contamination estimates.
 *
 *   radGenMcTrees / radRecMcTrees
 *     Optional radiative DVCS MC trees. Missing radiative files are skipped
 *     with warnings because the standard production workflow can use imported
 *     pass-1 radiative correction values instead of recalculating them.
 *
 * The extended overload additionally fills:
 *
 *   currentStudyGenMcTrees / currentStudyRecMcTrees
 *     All available current-dependent DVCS MC files found in:
 *
 *       /work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen
 *
 *     using the filename pattern consumed by dvcs_current_dependence.py:
 *
 *       gen_dvcsgen_<period>_<current>nA_<beam>MeV.root
 *       rec_dvcsgen_<period>_<current>nA_<beam>MeV.root
 *
 *     with "nobkg" interpreted as 0 nA.
 *
 *     The generated tags are:
 *
 *       DVCS_<PeriodTag>_<current>nA_gen_current
 *
 *     and the reconstructed tags are:
 *
 *       DVCS_<PeriodTag>_<current>nA_rec_current
 *
 *     where <PeriodTag> is one of:
 *
 *       Sp18_inb
 *       Sp18_out
 *       Fa18_inb
 *       Fa18_out
 *       Sp19_inb
 *
 * Return:
 *   true on success, false on any fatal load error.
 */
bool loadTrees(std::map<std::string, TTree*>& dataTrees,
               std::map<std::string, TTree*>& genMcTrees,
               std::map<std::string, TTree*>& recMcTrees,
               std::map<std::string, TTree*>& eppi0DataTrees,
               std::map<std::string, TTree*>& eppi0GenMcTrees,
               std::map<std::string, TTree*>& eppi0RecMcTrees,
               std::map<std::string, TTree*>& eppi0BkgTrees,
               std::map<std::string, TTree*>& radGenMcTrees,
               std::map<std::string, TTree*>& radRecMcTrees);

bool loadTrees(std::map<std::string, TTree*>& dataTrees,
               std::map<std::string, TTree*>& genMcTrees,
               std::map<std::string, TTree*>& recMcTrees,
               std::map<std::string, TTree*>& eppi0DataTrees,
               std::map<std::string, TTree*>& eppi0GenMcTrees,
               std::map<std::string, TTree*>& eppi0RecMcTrees,
               std::map<std::string, TTree*>& eppi0BkgTrees,
               std::map<std::string, TTree*>& radGenMcTrees,
               std::map<std::string, TTree*>& radRecMcTrees,
               std::map<std::string, TTree*>& currentStudyGenMcTrees,
               std::map<std::string, TTree*>& currentStudyRecMcTrees);

#endif // LOAD_TREES_H