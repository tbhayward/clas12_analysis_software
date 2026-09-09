#ifndef PI0_SYSTEMATICS_H
#define PI0_SYSTEMATICS_H

#include <string>

// Assign the pass-2 pi0-background subtraction systematic.
//
// The relative uncertainty on the acceptance-corrected pi0 background is
// inherited from the pass-1 dedicated variation study:
//     delta B_pi0 / B_pi0 = 0.072.
//
// The pass-2 contamination and acceptance-corrected yields are used to
// propagate that uncertainty to each final cross-section bin.  The routine
// also writes pass-1/pass-2 validation diagnostics.
bool pi0_systematics(
    const std::string& csv_path,
    const std::string& pass1_summary_path,
    const std::string& output_dir);

#endif
