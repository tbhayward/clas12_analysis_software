#include "charge_accumulation.h"
#include <set>
#include <vector>
#include "KinematicCuts.h"  // Assuming this is the correct path to your KinematicCuts class
#include "common_vars.h"    // Assuming this contains the global variables and necessary imports
#include "load_run_info_from_csv.h"

// This function accumulates the charges for the runs that pass the kinematic cuts.
void charge_accumulation(TTreeReader& dataReader, const std::vector<RunInfo>& run_info_list) {
    std::set<int> processedRuns; // To keep track of processed runs
    KinematicCuts kinematicCuts(dataReader); // Assumes KinematicCuts class is properly set up
    TTreeReaderValue<int> runnum(dataReader, "runnum"); // For retrieving the runnum from the data
    while (dataReader.Next()) {
        if (processedRuns.find(*runnum) == processedRuns.end()) {
            processedRuns.insert(*runnum);
            // Find run_info for the current run and update cpp, cpm, cmp, cmm
            for (const auto& run_info : run_info_list) {
                if (run_info.runnum == *runnum) {
                    if (run_info.target_polarization > 0) {
                        cpp += run_info.positive_charge;
                        cmp += run_info.negative_charge;
                    } else if (run_info.target_polarization < 0) {
                        cpm += run_info.positive_charge;
                        cmm += run_info.negative_charge;
                    }
                    break;
                }
            }
        }
    }
}
