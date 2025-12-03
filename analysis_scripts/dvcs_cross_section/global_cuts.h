#ifndef GLOBAL_CUTS_H
#define GLOBAL_CUTS_H

#include <string>
#include <vector>

/**
 * Global, universal event-level cuts for DVCS analysis.
 *
 * We apply these in ALL stages that iterate over events, before any
 * stage-specific logic:
 *    (-t1) < 1.0
 *    open_angle_ep2 > 5.0 deg
 *    pTmiss <= 0.20 (GeV)
 *
 * Keep this module simple and explicit so it is easy to read, copy, and reuse.
 */

struct GlobalCutConfig {
    // Require (-t1) < t1_abs_max. The analysis uses t1 with minus sign in mind,
    // so this expresses the absolute limit concisely.
    double t1_abs_max = 1.0;     // i.e. cut is: (-t1) < 1.0

    // Require open_angle_ep2 > open_angle_min_deg.
    double open_angle_min_deg = 5.0;

    // Require pTmiss <= pTmiss_max.
    double pTmiss_max = 0.20;    // GeV

    // NEW: runs to exclude globally (data quality blacklist).
    // Default-initialize with your known bad runs.
    std::vector<int> excluded_runs = {
        3867, 5046, 5047, 5051, 5128, 5129, 5130, 5160, 5163, 5165, 5166, 5167, 5168, 5169,
        5180, 5181, 5182, 5183, 5247, 5448, 5495, 5496, 5615, 5567
    };
};

// Returns a reference to a process-wide default configuration.
const GlobalCutConfig& default_global_cuts();

// Returns true if the triple (t1, open_angle_ep2_deg, pTmiss) passes the config.
bool passes_global_cuts(double t1,
                        double open_angle_ep2_deg,
                        double pTmiss,
                        const GlobalCutConfig& cfg = default_global_cuts());

// NEW: check if a run is globally blacklisted.
bool is_excluded_run(int runnum,
                     const GlobalCutConfig& cfg = default_global_cuts());

// Convenience: ROOT-style TCut string for Draw/Project workflows.
// Example: "(-t1) < 1.0 && open_angle_ep2 > 5.0 && pTmiss <= 0.20"
std::string global_cuts_tcut(const GlobalCutConfig& cfg = default_global_cuts());

// Persist the configuration once for provenance.
void write_global_cuts_config_json(const std::string& out_json_dir,
                                   const GlobalCutConfig& cfg = default_global_cuts());

#endif // GLOBAL_CUTS_H