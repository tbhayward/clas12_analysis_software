#include "kinematic_cuts.h"

// Function that applies the kinematic cuts
bool apply_kinematic_cuts(double t, double open_angle_ep2, double Emiss2, double Mx2_1, double pTmiss) {
    // Applying the specified cuts
    if (-t >= 1) return false;    
    if (open_angle_ep2 <= 10) return false; 
    if (Emiss2 >= 1) return false;        
    if (Mx2_1 >= 0.25) return false;        
    if (pTmiss >= 0.10) return false;       

    return true;  // Event passes all cuts
}