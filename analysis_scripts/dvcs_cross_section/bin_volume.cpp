#include "bin_volume.h"
#include <cmath>

double calculate_bin_volume(double xB_min, double xB_max,
                           double Q2_min, double Q2_max,
                           double t_min, double t_max, 
                           double phi_min, double phi_max,
                           double E_beam) {
    const int n_steps = 10;
    const double Mp = 0.938272; // Proton mass in GeV/c²
    int valid_count = 0;

    // Convert t to negative values (physical convention)
    double t_phys_min = -t_max;
    double t_phys_max = -t_min;

    // Sub-bin widths
    double dxB = (xB_max - xB_min)/n_steps;
    double dQ2 = (Q2_max - Q2_min)/n_steps;
    double dt = (t_phys_max - t_phys_min)/n_steps;

    for(int i=0; i<n_steps; i++) {
        double xB = xB_min + (i+0.5)*dxB;
        
        for(int j=0; j<n_steps; j++) {
            double Q2 = Q2_min + (j+0.5)*dQ2;
            
            // Calculate t_min for these kinematics
            // double sqrt_term = sqrt(1 + (4*Mp*Mp*xB*xB)/Q2);
            // double t_min_val = -Q2*(1-xB)*(1-xB)/(xB*(1 + sqrt_term));
            t_min_val = -pow((Mp*x),2)/(1-x);
            
            for(int k=0; k<n_steps; k++) {
                double t = t_phys_min + (k+0.5)*dt;
                
                // Calculate y and W
                double y = Q2/(2*Mp*xB*E_beam);
                double W = sqrt(Mp*Mp + Q2*(1.0/xB - 1));
                
                // Check kinematic cuts
                if(t > t_min_val && y > 0.19 && y < 0.8 && W > 2.0) {
                    valid_count++;
                }
            }
        }
    }

    double fraction = static_cast<double>(valid_count)/(n_steps*n_steps*n_steps);
    double geometric_volume = (xB_max - xB_min) * 
                             (Q2_max - Q2_min) * 
                             (t_max - t_min) * // Use original positive values
                             (phi_max - phi_min);
    
    return geometric_volume * fraction;
}