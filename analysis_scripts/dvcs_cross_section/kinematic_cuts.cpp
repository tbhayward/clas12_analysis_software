#include "kinematic_cuts.h"

// Function that applies the kinematic cuts
bool apply_kinematic_cuts(double t, double open_angle_ep2, double theta_neutral_neutral, double Emiss2, double Mx2_2, double pTmiss) {
    // Applying the specified cuts

    // dvcsgen cff inputs range
    if (-t >= 1) return false;  

    // remove radiative photons around e'
    if (open_angle_ep2 <= 5) return false; 

    // exclusivity cuts
    if (theta_neutral_neutral > 0.7) return false;
    if (Emiss2 >= 1) return false;        
    if (pTmiss >= 0.15) return false;
    if (Mx2_2 >= 1.75) return false;               

    return true;  // Event passes all cuts
}