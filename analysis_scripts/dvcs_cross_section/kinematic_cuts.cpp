#include "kinematic_cuts.h"

// Function that applies the kinematic cuts
bool apply_kinematic_cuts(double t, double open_angle_ep2, double Emiss2, double Mx2_1, double pTmiss) {
    // Applying the specified cuts
    if (t >= 1) return false;               // -t < 1
    if (open_angle_ep2 <= 10) return false; // open_angle_ep2 > 10
    if (Emiss2 >= 1) return false;          // Emiss2 < 1
    if (Mx2_1 >= 0.5) return false;         // Mx2_1 < 0.5
    if (pTmiss >= 0.2) return false;        // pTmiss < 0.2

    return true;  // Event passes all cuts
}