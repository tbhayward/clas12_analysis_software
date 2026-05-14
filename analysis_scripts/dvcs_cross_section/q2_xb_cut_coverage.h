#ifndef Q2_XB_CUT_COVERAGE_H
#define Q2_XB_CUT_COVERAGE_H

#include <map>
#include <string>

class TTree;

bool plot_q2_xb_cut_coverage(
    const std::map<std::string, TTree*>& dvcsDataTrees,
    const std::string& pass2_csv_path,
    const std::string& combined_cuts_json_path,
    const std::string& output_dir
);

#endif // Q2_XB_CUT_COVERAGE_H