#ifndef UNFOLDING_H
#define UNFOLDING_H

#include <string>

// Compute acceptance-corrected (unfolded) yields:
//   - Reads signal-yield triples and acceptance triples from the CSV
//   - Fills "acceptance corrected yield, ep->epg, exp, <label>, <hel>"
//     for periods and combined groups
//   - Produces xB-sliced canvases vs phi for pos/neg helicities,
//     saved under output/unfolded_yields/<PeriodOrGroupDir>/.
bool update_unfolded_yields_csv(const std::string& csv_path,
                                const std::string& out_root_dir);

#endif